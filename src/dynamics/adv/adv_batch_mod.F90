! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   The tracer advection can be seperated into different batches. Each batch
!   can have different time step size. The wind and mass flux are accumulated
!   along model integration, and averaged to middle time level of advection time
!   step cycle.
!
!   The batch type allocates necessary arrays, and provides wind accumulation
!   subroutines.
!
! Note:
!
!   - It needs to verify the wind and mass flux accumulation manners:
!     Averaging wind and mass flux on n + 1/2 time level, or n time level.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module adv_batch_mod

  use flogger
  use const_mod
  use namelist_mod
  use time_mod
  use allocator_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_parallel_mod
  use vert_coord_mod

  implicit none

  private

  public adv_batch_type

  ! Different tracers can be combined into one batch, and advected in different
  ! frequencies.
  type adv_batch_type
    character(30) :: scheme = 'N/A'
    character(10) :: loc  = 'cell'
    character(30) :: name = ''
    logical  :: dynamic   = .false.
    integer  :: ntracers  = 1
    integer  :: nstep     = 0       ! Number of dynamic steps for one adv step
    integer  :: step      = 0       ! Step counter
    real(r8) :: dt                  ! Advection time step size in seconds
    integer , allocatable :: idx(:) ! Global index of tracers in this batch
    type(latlon_field3d_type) old_m
    type(latlon_field3d_type) mfx
    type(latlon_field3d_type) mfy
    type(latlon_field3d_type) mz
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) we
    type(latlon_field3d_type) mfx0
    type(latlon_field3d_type) mfy0
    type(latlon_field3d_type) mx0
    type(latlon_field3d_type) my0
    type(latlon_field3d_type) mz0
    type(latlon_field3d_type) dmf
    type(latlon_field2d_type) dmgs
    type(latlon_field3d_type) qmf_lon
    type(latlon_field3d_type) qmf_lat
    type(latlon_field3d_type) qmf_lev
    ! FFSL variables
    type(latlon_field3d_type) cflx
    type(latlon_field3d_type) cfly
    type(latlon_field3d_type) cflz
    type(latlon_field3d_type) divx
    type(latlon_field3d_type) divy
    type(latlon_field3d_type) qx
    type(latlon_field3d_type) qy
  contains
    procedure :: init          => adv_batch_init
    procedure :: clear         => adv_batch_clear
    procedure :: copy_old_m    => adv_batch_copy_old_m
    procedure :: set_wind      => adv_batch_set_wind
    procedure :: accum_wind    => adv_batch_accum_wind
    procedure :: prepare       => adv_batch_prepare
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, filter_mesh, filter_halo, mesh, halo, scheme, batch_loc, batch_name, dt, dynamic, idx)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: scheme
    character(*), intent(in) :: batch_loc
    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    integer, intent(in), optional :: idx(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    this%scheme   = scheme
    this%loc      = batch_loc
    this%name     = batch_name
    this%dt       = dt
    this%dynamic  = dynamic
    this%nstep    = dt / dt_dyn
    this%step     = 0

    select case (batch_loc)
    case ('cell')
      if (.not. this%dynamic) then
        name      = trim(this%name) // '_old_m'
        long_name = 'Saved dry-air weight'
        units     = 'Pa'
        call this%old_m%init(batch_name, long_name, units, 'cell', mesh, halo)

        name      = trim(this%name) // '_mfx'
        long_name = 'Mass flux in x direction'
        units     = 'Pa m s-1'
        call this%mfx%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_mfy'
        long_name = 'Mass flux in y direction'
        units     = 'Pa m s-1'
        call this%mfy%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_mz'
        long_name = 'Dry-air weight on half level'
        units     = 'Pa'
        call this%mz%init(batch_name, long_name, units, 'lev', mesh, halo)

        name      = trim(this%name) // '_u'
        long_name = 'U wind component for advection of ' // trim(this%name)
        units     = 'm s-1'
        call this%u%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_v'
        long_name = 'V wind component for advection of ' // trim(this%name)
        units     = 'm s-1'
        call this%v%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_we'
        long_name = 'Vertical mass flux for advection of ' // trim(this%name)
        units     = 'Pa s-1'
        call this%we%init(batch_name, long_name, units, 'lev', mesh, halo)

        name      = trim(this%name) // '_mfx0'
        long_name = 'Mass flux in x direction'
        units     = 'Pa m s-1'
        call this%mfx0%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_mfy0'
        long_name = 'Mass flux in y direction'
        units     = 'Pa m s-1'
        call this%mfy0%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_mx0'
        long_name = ''
        units     = ''
        call this%mx0%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_my0'
        long_name = ''
        units     = ''
        call this%my0%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_mz0'
        long_name = ''
        units     = ''
        call this%mz0%init(batch_name, long_name, units, 'lev', mesh, halo)

        name      = trim(this%name) // '_dmf'
        long_name = ''
        units     = ''
        call this%dmf%init(batch_name, long_name, units, 'cell', filter_mesh, filter_halo)

        name      = trim(this%name) // '_dmgs'
        long_name = ''
        units     = ''
        call this%dmgs%init(batch_name, long_name, units, 'cell', mesh, halo)
      end if
      name        = trim(this%name) // '_qmf_lon'
      long_name   = 'Tracer mass flux in x direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lon%init(batch_name, long_name, units, 'lon', mesh, halo)

      name        = trim(this%name) // '_qmf_lat'
      long_name   = 'Tracer mass flux in y direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lat%init(batch_name, long_name, units, 'lat', mesh, halo)

      name        = trim(this%name) // '_qmf_lev'
      long_name   = 'Tracer mass flux in z direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lev%init(batch_name, long_name, units, 'lev', mesh, halo)
      select case (this%scheme)
      case ('ffsl')
        name      = trim(this%name) // '_cflx'
        long_name = 'CFL number in x direction'
        units     = ''
        call this%cflx%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_cfly'
        long_name = 'CFL number in y direction'
        units     = ''
        call this%cfly%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_cflz'
        long_name = 'CFL number in z direction'
        units     = ''
        call this%cflz%init(batch_name, long_name, units, 'lev', mesh, halo)

        name      = trim(this%name) // '_divx'
        long_name = 'Horizontal mass flux divergence in x direction'
        units     = 's-1'
        call this%divx%init(batch_name, long_name, units, 'lon', mesh, halo)

        name      = trim(this%name) // '_divy'
        long_name = 'Horizontal mass flux divergence in y direction'
        units     = 's-1'
        call this%divy%init(batch_name, long_name, units, 'lat', mesh, halo)

        name      = trim(this%name) // '_qx'
        long_name = 'Tracer mass after advection in x direction'
        units     = 'kg kg-1'
        call this%qx%init(batch_name, long_name, units, 'cell', filter_mesh, filter_halo, halo_cross_pole=.true.)

        name      = trim(this%name) // '_qy'
        long_name = 'Tracer mass after advection in y direction'
        units     = 'kg kg-1'
        call this%qy%init(batch_name, long_name, units, 'cell', filter_mesh, filter_halo)
      end select
    case ('lev')
      ! Only for nonhydrostatic dynamic calculation.
      name        = trim(this%name) // '_qmf_lon'
      long_name   = 'Tracer mass flux in x direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lon%init(batch_name, long_name, units, 'lev_lon', mesh, halo)

      name        = trim(this%name) // '_qmf_lat'
      long_name   = 'Tracer mass flux in y direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lat%init(batch_name, long_name, units, 'lev_lat', mesh, halo)

      name        = trim(this%name) // '_qmf_lev'
      long_name   = 'Tracer mass flux in z direction'
      units       = 'Pa kg kg-1 m s-1'
      call this%qmf_lev%init(batch_name, long_name, units, 'cell', mesh, halo)
      select case (this%scheme)
      case ('ffsl')
        name      = trim(this%name) // '_cflx'
        long_name = 'CFL number in x direction'
        units     = ''
        call this%cflx%init(batch_name, long_name, units, 'lev_lon', mesh, halo)

        name      = trim(this%name) // '_cfly'
        long_name = 'CFL number in y direction'
        units     = ''
        call this%cfly%init(batch_name, long_name, units, 'lev_lat', mesh, halo)

        name      = trim(this%name) // '_cflz'
        long_name = 'CFL number in z direction'
        units     = ''
        call this%cflz%init(batch_name, long_name, units, 'cell', mesh, halo)

        name      = trim(this%name) // '_divx'
        long_name = 'Horizontal mass flux divergence in x direction'
        units     = 's-1'
        call this%divx%init(batch_name, long_name, units, 'lev_lon', mesh, halo)

        name      = trim(this%name) // '_divy'
        long_name = 'Horizontal mass flux divergence in y direction'
        units     = 's-1'
        call this%divy%init(batch_name, long_name, units, 'lev_lat', mesh, halo)

        name      = trim(this%name) // '_qx'
        long_name = 'Tracer mass after advection in x direction'
        units     = 'kg kg-1'
        call this%qx%init(batch_name, long_name, units, 'lev', filter_mesh, filter_halo, halo_cross_pole=.true.)

        name      = trim(this%name) // '_qy'
        long_name = 'Tracer mass after advection in y direction'
        units     = 'kg kg-1'
        call this%qy%init(batch_name, long_name, units, 'lev', filter_mesh, filter_halo)
      end select
    case default
      call log_error('Invalid grid location ' // trim(batch_loc) // '!', __FILE__, __LINE__)
    end select

    if (present(idx)) then
      this%ntracers = size(idx)
      allocate(this%idx(this%ntracers))
      this%idx = idx
    end if

    call time_add_alert(batch_name, seconds=dt)

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    if (allocated (this%idx)) deallocate(this%idx)

    call this%old_m  %clear()
    call this%mfx    %clear()
    call this%mfy    %clear()
    call this%mz     %clear()
    call this%u      %clear()
    call this%v      %clear()
    call this%we     %clear()
    call this%mfx0   %clear()
    call this%mfy0   %clear()
    call this%mx0    %clear()
    call this%my0    %clear()
    call this%mz0    %clear()
    call this%dmf    %clear()
    call this%dmgs   %clear()
    call this%qmf_lon%clear()
    call this%qmf_lat%clear()
    call this%qmf_lev%clear()
    call this%cflx   %clear()
    call this%cfly   %clear()
    call this%cflz   %clear()
    call this%divx   %clear()
    call this%divy   %clear()
    call this%qx     %clear()
    call this%qy     %clear()

    this%loc         = 'cell'
    this%name        = ''
    this%dt          = 0
    this%dynamic     = .false.
    this%ntracers    = 0
    this%nstep       = 0
    this%step        = 0

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_old_m(this, m)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: m

    this%old_m%d = m%d

  end subroutine adv_batch_copy_old_m

  subroutine adv_batch_set_wind(this, u_lon, v_lat, we_lev, mfx_lon, mfy_lat, dmg_lev)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u_lon
    type(latlon_field3d_type), intent(in) :: v_lat
    type(latlon_field3d_type), intent(in) :: we_lev
    type(latlon_field3d_type), intent(in) :: mfx_lon
    type(latlon_field3d_type), intent(in) :: mfy_lat
    type(latlon_field3d_type), intent(in) :: dmg_lev

    call this%u  %link(u_lon  )
    call this%v  %link(v_lat  )
    call this%we %link(we_lev )
    call this%mfx%link(mfx_lon)
    call this%mfy%link(mfy_lat)
    call this%mz %link(dmg_lev)

    if (this%scheme == 'ffsl') call this%prepare()

  end subroutine adv_batch_set_wind

  subroutine adv_batch_accum_wind(this, dmg_lon, dmg_lat, dmg_lev, mfx_lon, mfy_lat)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: dmg_lon
    type(latlon_field3d_type), intent(in) :: dmg_lat
    type(latlon_field3d_type), intent(in) :: dmg_lev
    type(latlon_field3d_type), intent(in) :: mfx_lon
    type(latlon_field3d_type), intent(in) :: mfy_lat

    real(r8) work(this%u%mesh%full_ids:this%u%mesh%full_ide,this%u%mesh%full_nlev)
    real(r8) pole(this%u%mesh%full_nlev)
    integer i, j, k

    if (this%step == -1) then
      ! Reset step.
      this%mfx%d = this%mfx0%d
      this%mfy%d = this%mfy0%d
      this%u  %d = this%mx0 %d
      this%v  %d = this%my0 %d
      this%mz %d = this%mz0 %d
      this%step = 1
    end if
    if (this%step == 0) then
      ! This is the first step.
      this%mfx%d = mfx_lon%d
      this%mfy%d = mfy_lat%d
      this%u  %d = dmg_lon%d
      this%v  %d = dmg_lat%d
      this%mz %d = dmg_lev%d
    else if (this%step == this%nstep) then
      ! This is the end step.
      this%mfx %d = (this%mfx%d + mfx_lon%d) / (this%nstep + 1)
      this%mfy %d = (this%mfy%d + mfy_lat%d) / (this%nstep + 1)
      this%u   %d = (this%u  %d + dmg_lon%d) / (this%nstep + 1)
      this%v   %d = (this%v  %d + dmg_lat%d) / (this%nstep + 1)
      this%mz  %d = (this%mz %d + dmg_lev%d) / (this%nstep + 1)
      this%mfx0%d = mfx_lon%d
      this%mfy0%d = mfy_lat%d
      this%mx0 %d = dmg_lon%d
      this%my0 %d = dmg_lat%d
      this%mz0 %d = dmg_lev%d
    else
      ! Accumulating.
      this%mfx %d = this%mfx%d + mfx_lon%d
      this%mfy %d = this%mfy%d + mfy_lat%d
      this%u   %d = this%u  %d + dmg_lon%d
      this%v   %d = this%v  %d + dmg_lat%d
      this%mz  %d = this%mz %d + dmg_lev%d
    end if
    this%step = merge(0, this%step + 1, this%dynamic)
    if (this%dynamic .or. this%step > this%nstep) then
      if (.not. this%dynamic) this%step = -1
      associate (mesh => this%u%mesh, &
                 dt   => this%dt    , &
                 mfx  => this%mfx   , &
                 mfy  => this%mfy   , &
                 mx   => this%u     , &
                 my   => this%v     , &
                 mz   => this%mz    , &
                 u    => this%u     , &
                 v    => this%v     , &
                 we   => this%we    , &
                 dmf  => this%dmf   , &
                 dmgs => this%dmgs  )
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            u%d(i,j,k) = mfx%d(i,j,k) / mx%d(i,j,k)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = mfy%d(i,j,k) / my%d(i,j,k)
          end do
        end do
      end do
      ! Diagnose horizontal mass flux divergence.
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            dmf%d(i,j,k) = ((                     &
              mfx%d(i,j,k) - mfx%d(i-1,j,k)       &
            ) * mesh%le_lon(j) + (                &
              mfy%d(i,j  ,k) * mesh%le_lat(j  ) - &
              mfy%d(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = mfy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            dmf%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = -mfy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            dmf%d(i,j,k) = pole(k)
          end do
        end do
      end if
      ! Diagnose surface dry air pressure tendency.
      dmgs%d = 0
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dmgs%d(i,j) = dmgs%d(i,j) - dmf%d(i,j,k)
          end do
        end do
      end do
      ! Diagnose vertical mass flux.
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            we%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum(dmf%d(i,j,1:k-1))
          end do
        end do
      end do
      end associate
      if (this%scheme == 'ffsl') call this%prepare()
    end if

  end subroutine adv_batch_accum_wind

  subroutine adv_batch_prepare(this)

    class(adv_batch_type), intent(inout) :: this

    real(r8) work(this%u%mesh%full_ids:this%u%mesh%full_ide,this%u%mesh%half_nlev)
    real(r8) pole(this%u%mesh%half_nlev)
    integer ks, ke, i, j, k

    associate (mesh => this%u%mesh, &
               dt   => this%dt    , &
               mfx  => this%mfx   , &
               mfy  => this%mfy   , &
               mz   => this%mz    , &
               u    => this%u     , &
               v    => this%v     , &
               we   => this%we    , &
               cflx => this%cflx  , &
               cfly => this%cfly  , &
               cflz => this%cflz  , &
               divx => this%divx  , &
               divy => this%divy  )
    select case (this%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, this%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, this%loc == 'cell')
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            cflx%d(i,j,k) = u%d(i,j,k) * dt / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            cfly%d(i,j,k) = v%d(i,j,k) * dt / mesh%de_lat(j)
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            divx%d(i,j,k) = (u%d(i,j,k) - u%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
            divy%d(i,j,k) = (v%d(i,j,k) * mesh%le_lat(j) - v%d(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = v%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = -v%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('vtx')
    end select

    select case (this%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            cflz%d(i,j,k) = we%d(i,j,k) / mz%d(i,j,k) * dt
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            cflz%d(i,j,k) = we%d(i,j,k) / mz%d(i,j,k) * dt
          end do
        end do
      end do
    end select
    end associate

  end subroutine adv_batch_prepare

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
