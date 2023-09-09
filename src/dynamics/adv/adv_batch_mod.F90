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
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use latlon_parallel_mod
  use vert_coord_mod

  implicit none

  private

  public adv_batch_type

  ! Different tracers can be combined into one batch, and advected in different
  ! frequencies.
  type adv_batch_type
    type(mesh_type), pointer :: filter_mesh => null()
    type(mesh_type), pointer :: mesh => null()
    character(10) :: loc  = 'cell'
    character(30) :: name = ''
    logical  :: dynamic   = .false.
    integer  :: ntracers  = 1
    integer  :: nstep     = 0 ! Number of dynamic steps for one adv step
    integer  :: step      = 0 ! Step counter
    real(r8) :: dt            ! Advection time step size in seconds
    integer , allocatable, dimension(:    ) :: idx   ! Global index of tracers in this batch
    real(r8), allocatable, dimension(:,:,:) :: old_m ! Recorded old mass for converting mixing ratio
    real(r8), pointer    , dimension(:,:,:) :: mfx => null()
    real(r8), pointer    , dimension(:,:,:) :: mfy => null()
    real(r8), pointer    , dimension(:,:,:) :: mz  => null()
    real(r8), pointer    , dimension(:,:,:) :: u   => null()
    real(r8), pointer    , dimension(:,:,:) :: v   => null()
    real(r8), pointer    , dimension(:,:,:) :: we  => null()
    real(r8), allocatable, dimension(:,:,:) :: mfx0
    real(r8), allocatable, dimension(:,:,:) :: mfy0
    real(r8), allocatable, dimension(:,:,:) :: mx0
    real(r8), allocatable, dimension(:,:,:) :: my0
    real(r8), allocatable, dimension(:,:,:) :: mz0
    real(r8), allocatable, dimension(:,:,:) :: dmf
    real(r8), allocatable, dimension(:,:  ) :: dmgs
    ! The following arrays could be reused by different batches.
    real(r8), allocatable, dimension(:,:,:) :: qmf_lon
    real(r8), allocatable, dimension(:,:,:) :: qmf_lat
    real(r8), allocatable, dimension(:,:,:) :: qmf_lev
    ! FFSL variables
    real(r8), allocatable, dimension(:,:,:) :: cflx ! CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cfly ! CFL number along y-axis
    real(r8), allocatable, dimension(:,:,:) :: cflz ! CFL number along z-axis
    real(r8), allocatable, dimension(:,:,:) :: divx ! Divergence along x-axis
    real(r8), allocatable, dimension(:,:,:) :: divy ! Divergence along y-axis
    real(r8), allocatable, dimension(:,:,:) :: qx   ! Tracer mixing ratio due to advective operator along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy   ! Tracer mixing ratio due to advective operator along y axis
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

  subroutine adv_batch_init(this, filter_mesh, mesh, loc, name, dt, dynamic, idx)

    class(adv_batch_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: filter_mesh
    type(mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: loc
    character(*), intent(in) :: name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    integer, intent(in), optional :: idx(:)

    call this%clear()

    this%filter_mesh => filter_mesh
    this%mesh        => mesh
    this%loc         = loc
    this%name        = name
    this%dt          = dt
    this%dynamic     = dynamic
    this%nstep       = dt / dt_dyn
    this%step        = 0

    select case (loc)
    case ('cell')
      if (.not. this%dynamic) then
        call allocate_array(mesh, this%old_m, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mfx  , half_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mfy  , full_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mz   , full_lon=.true., full_lat=.true., half_lev=.true.)
        call allocate_array(mesh, this%u    , half_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%v    , full_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%we   , full_lon=.true., full_lat=.true., half_lev=.true.)
        call allocate_array(mesh, this%mfx0 , half_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mfy0 , full_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mx0  , half_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%my0  , full_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%mz0  , full_lon=.true., full_lat=.true., half_lev=.true.)
        call allocate_array(mesh, this%dmf  , full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%dmgs , full_lon=.true., full_lat=.true.)
      end if
      call allocate_array(mesh, this%qmf_lon, half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qmf_lat, full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qmf_lev, full_lon=.true., full_lat=.true., half_lev=.true.)
      select case (adv_scheme)
      case ('ffsl')
        call allocate_array(mesh, this%cflx , half_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%cfly , full_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%cflz , full_lon=.true., full_lat=.true., half_lev=.true.)
        call allocate_array(mesh, this%divx , full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%divy , full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(filter_mesh, this%qx, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(filter_mesh, this%qy, full_lon=.true., full_lat=.true., full_lev=.true.)
      end select
    case default
      call log_error('Invalid grid location ' // trim(loc) // '!', __FILE__, __LINE__)
    end select

    if (present(idx)) then
      this%ntracers = size(idx)
      allocate(this%idx(this%ntracers))
      this%idx = idx
    end if

    call time_add_alert(name, seconds=dt/time_scale)

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    if (allocated (this%idx    )) deallocate(this%idx    )
    if (allocated (this%old_m  )) deallocate(this%old_m  )
    if (.not. this%dynamic) then
      if (associated(this%mfx  )) deallocate(this%mfx    )
      if (associated(this%mfy  )) deallocate(this%mfy    )
      if (associated(this%mz   )) deallocate(this%mz     )
      if (associated(this%u    )) deallocate(this%u      )
      if (associated(this%v    )) deallocate(this%v      )
      if (associated(this%we   )) deallocate(this%we     )
    else
      this%mfx => null()
      this%mfy => null()
      this%mz  => null()
      this%u   => null()
      this%v   => null()
      this%we  => null()
    end if
    if (allocated (this%mfx0   )) deallocate(this%mfx0   )
    if (allocated (this%mfy0   )) deallocate(this%mfy0   )
    if (allocated (this%mx0    )) deallocate(this%mx0    )
    if (allocated (this%my0    )) deallocate(this%my0    )
    if (allocated (this%mz0    )) deallocate(this%mz0    )
    if (allocated (this%dmf    )) deallocate(this%dmf    )
    if (allocated (this%dmgs   )) deallocate(this%dmgs   )
    if (allocated (this%qmf_lon)) deallocate(this%qmf_lon)
    if (allocated (this%qmf_lat)) deallocate(this%qmf_lat)
    if (allocated (this%qmf_lev)) deallocate(this%qmf_lev)
    if (allocated (this%cflx   )) deallocate(this%cflx   )
    if (allocated (this%cfly   )) deallocate(this%cfly   )
    if (allocated (this%cflz   )) deallocate(this%cflz   )
    if (allocated (this%divx   )) deallocate(this%divx   )
    if (allocated (this%divy   )) deallocate(this%divy   )
    if (allocated (this%qx     )) deallocate(this%qx     )
    if (allocated (this%qy     )) deallocate(this%qy     )

    this%filter_mesh => null()
    this%mesh        => null()
    this%loc         = 'cell'
    ! this%name        = ''
    this%dt          = 0
    this%dynamic     = .false.
    this%ntracers    = 0
    this%nstep       = 0
    this%step        = 0

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_old_m(this, m)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: m(:,:,:) ! Assume caller knows what the array shape is.

    this%old_m = m

  end subroutine adv_batch_copy_old_m

  subroutine adv_batch_set_wind(this, u_lon, v_lat, we_lev, mfx_lon, mfy_lat, dmg_lev)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in), target, contiguous :: u_lon  (:,:,:)
    real(r8), intent(in), target, contiguous :: v_lat  (:,:,:)
    real(r8), intent(in), target, contiguous :: we_lev (:,:,:)
    real(r8), intent(in), target, contiguous :: mfx_lon(:,:,:)
    real(r8), intent(in), target, contiguous :: mfy_lat(:,:,:)
    real(r8), intent(in), target, contiguous :: dmg_lev(:,:,:)

    this%u  (this%mesh%half_ims:this%mesh%half_ime, &
             this%mesh%full_jms:this%mesh%full_jme, &
             this%mesh%full_kms:this%mesh%full_kme) => u_lon
    this%v  (this%mesh%full_ims:this%mesh%full_ime, &
             this%mesh%half_jms:this%mesh%half_jme, &
             this%mesh%full_kms:this%mesh%full_kme) => v_lat
    this%we (this%mesh%full_ims:this%mesh%full_ime, &
             this%mesh%full_jms:this%mesh%full_jme, &
             this%mesh%half_kms:this%mesh%half_kme) => we_lev
    this%mfx(this%mesh%half_ims:this%mesh%half_ime, &
             this%mesh%full_jms:this%mesh%full_jme, &
             this%mesh%full_kms:this%mesh%full_kme) => mfx_lon
    this%mfy(this%mesh%full_ims:this%mesh%full_ime, &
             this%mesh%half_jms:this%mesh%half_jme, &
             this%mesh%full_kms:this%mesh%full_kme) => mfy_lat
    this%mz (this%mesh%full_ims:this%mesh%full_ime, &
             this%mesh%full_jms:this%mesh%full_jme, &
             this%mesh%half_kms:this%mesh%half_kme) => dmg_lev

    call this%prepare()

  end subroutine adv_batch_set_wind

  subroutine adv_batch_accum_wind(this, dmg_lon, dmg_lat, dmg_lev, mfx_lon, mfy_lat)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dmg_lon(this%mesh%half_ims:this%mesh%half_ime, &
                                    this%mesh%full_jms:this%mesh%full_jme, &
                                    this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in) :: dmg_lat(this%mesh%full_ims:this%mesh%full_ime, &
                                    this%mesh%half_jms:this%mesh%half_jme, &
                                    this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in) :: dmg_lev(this%mesh%full_ims:this%mesh%full_ime, &
                                    this%mesh%full_jms:this%mesh%full_jme, &
                                    this%mesh%half_kms:this%mesh%half_kme)
    real(r8), intent(in) :: mfx_lon(this%mesh%half_ims:this%mesh%half_ime, &
                                    this%mesh%full_jms:this%mesh%full_jme, &
                                    this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in) :: mfy_lat(this%mesh%full_ims:this%mesh%full_ime, &
                                    this%mesh%half_jms:this%mesh%half_jme, &
                                    this%mesh%full_kms:this%mesh%full_kme)

    real(r8) work(this%mesh%full_ids:this%mesh%full_ide,this%mesh%full_nlev)
    real(r8) pole(this%mesh%full_nlev)
    integer i, j, k

    if (this%step == -1) then
      ! Reset step.
      this%mfx = this%mfx0
      this%mfy = this%mfy0
      this%u   = this%mx0
      this%v   = this%my0
      this%mz  = this%mz0
      this%step = 1
    end if
    if (this%step == 0) then
      ! This is the first step.
      this%mfx = mfx_lon
      this%mfy = mfy_lat
      this%u   = dmg_lon
      this%v   = dmg_lat
      this%mz  = dmg_lev
    else if (this%step == this%nstep) then
      ! This is the end step.
      this%mfx = (this%mfx + mfx_lon) / (this%nstep + 1)
      this%mfy = (this%mfy + mfy_lat) / (this%nstep + 1)
      this%u   = (this%u   + dmg_lon) / (this%nstep + 1)
      this%v   = (this%v   + dmg_lat) / (this%nstep + 1)
      this%mz  = (this%mz  + dmg_lev) / (this%nstep + 1)
      this%mfx0 = mfx_lon
      this%mfy0 = mfy_lat
      this%mx0  = dmg_lon
      this%my0  = dmg_lat
      this%mz0  = dmg_lev
    else
      ! Accumulating.
      this%mfx = this%mfx + mfx_lon
      this%mfy = this%mfy + mfy_lat
      this%u   = this%u   + dmg_lon
      this%v   = this%v   + dmg_lat
      this%mz  = this%mz  + dmg_lev
    end if
    this%step = merge(0, this%step + 1, this%dynamic)
    if (this%dynamic .or. this%step > this%nstep) then
      if (.not. this%dynamic) this%step = -1
      associate (mesh => this%mesh, &
                 dt   => this%dt  , &
                 mfx  => this%mfx , &
                 mfy  => this%mfy , &
                 mx   => this%u   , &
                 my   => this%v   , &
                 mz   => this%mz  , &
                 u    => this%u   , &
                 v    => this%v   , &
                 we   => this%we  , &
                 dmf  => this%dmf , &
                 dmgs => this%dmgs)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            u(i,j,k) = mfx(i,j,k) / mx(i,j,k)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v(i,j,k) = mfy(i,j,k) / my(i,j,k)
          end do
        end do
      end do
      ! Diagnose horizontal mass flux divergence.
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            dmf(i,j,k) = ((                     &
              mfx(i,j,k) - mfx(i-1,j,k)         &
            ) * mesh%le_lon(j) + (              &
              mfy(i,j  ,k) * mesh%le_lat(j  ) - &
              mfy(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = mfy(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            dmf(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = -mfy(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            dmf(i,j,k) = pole(k)
          end do
        end do
      end if
      ! Diagnose surface dry air pressure tendency.
      dmgs = 0
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dmgs(i,j) = dmgs(i,j) - dmf(i,j,k)
          end do
        end do
      end do
      ! Diagnose vertical mass flux.
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            we(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs(i,j)) - sum(dmf(i,j,1:k-1))
          end do
        end do
      end do
      end associate
      call this%prepare()
    end if

  end subroutine adv_batch_accum_wind

  subroutine adv_batch_prepare(this)

    class(adv_batch_type), intent(inout) :: this

    real(r8) work(this%mesh%full_ids:this%mesh%full_ide,this%mesh%full_nlev)
    real(r8) pole(this%mesh%full_nlev)
    integer i, j, k

    associate (mesh => this%mesh, &
               dt   => this%dt  , &
               mfx  => this%mfx , &
               mfy  => this%mfy , &
               mz   => this%mz  , &
               u    => this%u   , &
               v    => this%v   , &
               we   => this%we  , &
               cflx => this%cflx, &
               cfly => this%cfly, &
               cflz => this%cflz, &
               divx => this%divx, &
               divy => this%divy)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          cflx(i,j,k) = u(i,j,k) * dt / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          cfly(i,j,k) = v(i,j,k) * dt / mesh%de_lat(j)
        end do
      end do
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          divx(i,j,k) = (u(i,j,k) - u(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
          divy(i,j,k) = (v(i,j,k) * mesh%le_lat(j) - v(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          divy(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -v(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          divy(i,j,k) = pole(k)
        end do
      end do
    end if
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          cflz(i,j,k) = we(i,j,k) / mz(i,j,k) * dt
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_prepare

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
