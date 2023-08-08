module adv_batch_mod

  use flogger
  use const_mod
  use namelist_mod
  use latlon_mesh_mod
  use time_mod
  use allocator_mod
  use latlon_parallel_mod

  implicit none

  private

  public adv_batch_type

  ! Different tracers can be combined into one batch, and adved in different frequencfly.
  type adv_batch_type
    type(latlon_mesh_type), pointer :: filter_mesh => null()
    type(latlon_mesh_type), pointer :: mesh => null()
    character(10) :: loc  = 'cell'
    character(30) :: name = ''
    logical  :: dynamic   = .false.
    integer  :: ntracers  = 1
    integer  :: nstep     = 0 ! Number of dynamic steps for one adv step
    integer  :: uv_step   = 0 ! Step counter for u and v
    integer  :: we_step   = 0 ! Step counter for we
    integer  :: mf_step   = 0 ! Step counter for mass flux
    real(r8) :: dt            ! Advection time step size in seconds
    integer , allocatable, dimension(:    ) :: idx   ! Index of tracers in this batch
    real(r8), allocatable, dimension(:,:,:) :: old_m ! Recorded old mass for converting mixing ratio
    real(r8), allocatable, dimension(:,:,:) :: mfx
    real(r8), allocatable, dimension(:,:,:) :: mfy
    real(r8), allocatable, dimension(:,:,:) :: m , m0
    real(r8), allocatable, dimension(:,:,:) :: u , u0
    real(r8), allocatable, dimension(:,:,:) :: v , v0
    real(r8), allocatable, dimension(:,:,:) :: we, we0
    real(r8), allocatable, dimension(:,:,:) :: we_imp
    real(r8), allocatable, dimension(:,:,:) :: cflx ! CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cfly ! CFL number along y-axis
    real(r8), allocatable, dimension(:,:,:) :: cflz ! CFL number along z-axis
    real(r8), allocatable, dimension(:,:,:) :: divx ! Divergence along x-axis
    real(r8), allocatable, dimension(:,:,:) :: divy ! Divergence along y-axis
    ! The following arrays could be reused by different batches.
    real(r8), allocatable, dimension(:,:,:) :: qmf_lon
    real(r8), allocatable, dimension(:,:,:) :: qmf_lat
    real(r8), allocatable, dimension(:,:,:) :: qmf_lev
    real(r8), allocatable, dimension(:,:,:) :: qx   ! Tracer mixing ratio due to advective operator along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy   ! Tracer mixing ratio due to advective operator along y axis
  contains
    procedure :: init          => adv_batch_init
    procedure :: clear         => adv_batch_clear
    procedure :: copy_old_m    => adv_batch_copy_old_m
    procedure :: accum_uv_cell => adv_batch_accum_uv_cell
    procedure :: accum_mf_cell => adv_batch_accum_mf_cell
    procedure :: accum_we_lev  => adv_batch_accum_we_lev
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, filter_mesh, mesh, loc, name, dt, dynamic, idx)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: loc
    character(*), intent(in) :: name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    integer, intent(in), optional :: idx(:)

    call this%clear()

    this%filter_mesh => filter_mesh
    this%mesh      => mesh
    this%loc       = loc
    this%name      = name
    this%dt        = dt
    this%dynamic   = dynamic
    this%nstep     = dt / dt_dyn
    this%uv_step   = 0
    this%we_step   = 0
    this%mf_step   = 0

    select case (loc)
    case ('cell')
      call allocate_array(mesh, this%old_m  , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfx    , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfy    , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%m      , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%m0     , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%u      , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%u0     , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%v      , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%v0     , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%we     , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%we0    , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%cflx   , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cfly   , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cflz   , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%divx   , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%divy   , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qmf_lon, half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qmf_lat, full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qmf_lev, full_lon=.true., full_lat=.true., half_lev=.true.)
      select case (adv_scheme)
      case ('ffsl')
        call allocate_array(filter_mesh, this%qx, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(filter_mesh, this%qy, full_lon=.true., full_lat=.true., full_lev=.true.)
      end select
    case default
      call log_error('Invalid grid location ' // trim(loc) // '!', __FILE__, __LINE__)
    end select
    if (use_ieva) then
      call allocate_array(mesh, this%we_imp , full_lon=.true., full_lat=.true., half_lev=.true.)
    end if

    if (present(idx)) then
      this%ntracers = size(idx)
      allocate(this%idx(this%ntracers))
      this%idx = idx
    end if

    call time_add_alert(name, seconds=dt/time_scale)

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    this%mesh      => null()
    this%loc       = 'cell'
    this%name      = ''
    this%dt        = 0
    this%dynamic   = .false.
    this%ntracers  = 0
    this%nstep     = 0
    this%uv_step   = 0
    this%we_step   = 0
    this%mf_step   = 0

    if (allocated(this%idx    )) deallocate(this%idx    )
    if (allocated(this%old_m  )) deallocate(this%old_m  )
    if (allocated(this%mfx    )) deallocate(this%mfx    )
    if (allocated(this%mfy    )) deallocate(this%mfy    )
    if (allocated(this%m      )) deallocate(this%m      )
    if (allocated(this%u      )) deallocate(this%u      )
    if (allocated(this%u0     )) deallocate(this%u0     )
    if (allocated(this%v      )) deallocate(this%v      )
    if (allocated(this%v0     )) deallocate(this%v0     )
    if (allocated(this%we     )) deallocate(this%we     )
    if (allocated(this%we0    )) deallocate(this%we0    )
    if (allocated(this%we_imp )) deallocate(this%we_imp )
    if (allocated(this%cflx   )) deallocate(this%cflx   )
    if (allocated(this%cfly   )) deallocate(this%cfly   )
    if (allocated(this%cflz   )) deallocate(this%cflz   )
    if (allocated(this%divx   )) deallocate(this%divx   )
    if (allocated(this%divy   )) deallocate(this%divy   )
    if (allocated(this%qmf_lon)) deallocate(this%qmf_lon)
    if (allocated(this%qmf_lat)) deallocate(this%qmf_lat)
    if (allocated(this%qmf_lev)) deallocate(this%qmf_lev)
    if (allocated(this%qx     )) deallocate(this%qx     )
    if (allocated(this%qy     )) deallocate(this%qy     )

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_old_m(this, m)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: m(this%mesh%full_ims:this%mesh%full_ime, &
                              this%mesh%full_jms:this%mesh%full_jme, &
                              this%mesh%full_kms:this%mesh%full_kme)

    this%old_m = m

  end subroutine adv_batch_copy_old_m

  subroutine adv_batch_accum_uv_cell(this, u, v, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: u(this%mesh%half_ims:this%mesh%half_ime, &
                              this%mesh%full_jms:this%mesh%full_jme, &
                              this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in) :: v(this%mesh%full_ims:this%mesh%full_ime, &
                              this%mesh%half_jms:this%mesh%half_jme, &
                              this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in), optional :: dt

    real(r8) work(this%mesh%full_ids:this%mesh%full_ide,this%mesh%full_nlev)
    real(r8) pole(this%mesh%full_nlev)
    real(r8) dt_
    real(r8) dx, x0, x1, x2, x3, u1, u2, u3, u4
    real(r8) dy, y0, y1, y2, y3, v1, v2, v3, v4
    integer i, j, k, l

    dt_ = merge(dt, this%dt, present(dt))

    associate (mesh => this%mesh)
    if (this%uv_step == -1) then
      this%u = this%u0
      this%v = this%v0
      this%uv_step = 1
    end if
    if (this%uv_step == 0) then
      this%u = u
      this%v = v
    else if (this%uv_step == this%nstep) then
      this%u = (this%u + u) / (this%nstep + 1)
      this%v = (this%v + v) / (this%nstep + 1)
      this%u0 = u
      this%v0 = v
    else
      this%u = this%u + u
      this%v = this%v + v
    end if
    this%uv_step = merge(0, this%uv_step + 1, this%dynamic)
    if (this%dynamic .or. this%uv_step > this%nstep) then
      if (.not. this%dynamic) this%uv_step = -1
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          dx = mesh%de_lon(j)
          do i = mesh%half_ids, mesh%half_ide
            this%cflx(i,j,k) = this%u(i,j,k) * dt_ / dx
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          dy = mesh%de_lat(j)
          do i = mesh%full_ids, mesh%full_ide
            this%cfly(i,j,k) = this%v(i,j,k) * dt_ / dy
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            this%divx(i,j,k) = (this%u(i,j,k) - this%u(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
            this%divy(i,j,k) = (this%v(i,j  ,k) * mesh%le_lat(j  ) - &
                                this%v(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = this%v(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            this%divy(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = -this%v(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%full_ids, mesh%full_ide
            this%divy(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    end associate

  end subroutine adv_batch_accum_uv_cell

  subroutine adv_batch_accum_mf_cell(this, mfx, mfy)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: mfx(this%mesh%half_ims:this%mesh%half_ime, &
                                this%mesh%full_jms:this%mesh%full_jme, &
                                this%mesh%full_kms:this%mesh%full_kme)
    real(r8), intent(in) :: mfy(this%mesh%full_ims:this%mesh%full_ime, &
                                this%mesh%half_jms:this%mesh%half_jme, &
                                this%mesh%full_kms:this%mesh%full_kme)

    if (this%mf_step == -1) then
      this%mfx = 0
      this%mfy = 0
      this%mf_step = 1
    end if
    if (this%mf_step == 0) then
      this%mfx = mfx
      this%mfy = mfy
    else if (this%mf_step == this%nstep) then
      this%mfx = (this%mfx + mfx) / this%nstep
      this%mfy = (this%mfy + mfy) / this%nstep
    else
      this%mfx = this%mfx + mfx
      this%mfy = this%mfy + mfy
    end if
    this%mf_step = merge(0, this%mf_step + 1, this%dynamic)
    if (.not. this%dynamic .and. this%mf_step > this%nstep) this%mf_step = -1

  end subroutine adv_batch_accum_mf_cell

  subroutine adv_batch_accum_we_lev(this, we, m, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: we(this%mesh%full_ims:this%mesh%full_ime, &
                               this%mesh%full_jms:this%mesh%full_jme, &
                               this%mesh%half_kms:this%mesh%half_kme)
    real(r8), intent(in) :: m (this%mesh%full_ims:this%mesh%full_ime, &
                               this%mesh%full_jms:this%mesh%full_jme, &
                               this%mesh%half_kms:this%mesh%half_kme)
    real(r8), intent(in), optional :: dt

    real(r8) dt_
    real(r8) z0, z1, z2, z3, w1, w2, w3, w4, deta
    integer i, j, k, l, ks, ke, s
    real(r8), parameter :: alpha_max = 1.1, &
                           alpha_min = 0.8, &
                           kesi      = 0.9
    real(r8) alpha_h, alpha_v, alpha_star_max, alpha_star_min, beta
    real(r8) work(this%mesh%full_ids:this%mesh%full_ide,this%mesh%full_nlev)
    real(r8) pole(this%mesh%full_nlev)

    dt_ = merge(dt, this%dt, present(dt))

    associate (mesh => this%mesh)
    if (this%we_step == -1) then
      this%we = this%we0
      this%m  = this%m0
      this%we_step = 1
    end if
    if (this%we_step == 0) then
      this%we = we
      this%m  = m
    else if (this%we_step == this%nstep) then
      this%we = (this%we + we) / (this%nstep + 1)
      this%m  = (this%m  + m ) / (this%nstep + 1)
      this%we0 = we
      this%m0  = m
    else
      this%we = this%we + we
      this%m  = this%m  + m
    end if
    this%we_step = merge(0, this%we_step + 1, this%dynamic)
    if (this%dynamic .or. this%we_step > this%nstep) then
      if (.not. this%dynamic) this%we_step = -1
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            this%cflz(i,j,k) = this%we(i,j,k) / this%m(i,j,k) * dt_
          end do
        end do
      end do
      if (use_ieva) then
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              if (this%we(i,j,k) >= 0.0) then
                alpha_h = dt_ * ((max(this%u(i,j  ,k-1), 0.) - min(this%u(i-1,j,k-1), 0.)) * mesh%le_lon(j) + &
                                 (max(this%v(i,j  ,k-1), 0.) * mesh%le_lat(j  ) - &
                                  min(this%v(i,j-1,k-1), 0.) * mesh%le_lat(j-1))) / mesh%area_cell(j)
              else
                alpha_h = dt_ * ((max(this%u(i,j  ,k), 0.) - min(this%u(i-1,j,k), 0.)) * mesh%le_lon(j) + &
                                 (max(this%v(i,j  ,k), 0.) * mesh%le_lat(j  ) - &
                                  min(this%v(i,j-1,k), 0.) * mesh%le_lat(j-1))) / mesh%area_cell(j)
              end if
              alpha_star_max = alpha_max - kesi * alpha_h
              alpha_star_min = alpha_min * alpha_star_max / alpha_max
              alpha_v = dt_ * abs(this%we(i,j,k) / this%m(i,j,k))
              if (alpha_v <= alpha_star_min) then
                beta = 1
              else if (alpha_v > alpha_star_min .and. alpha_v <= 2 * alpha_star_max - alpha_star_min) then
                beta = 1 / (1 + (alpha_v - alpha_star_min)**2 / (4 * alpha_star_max * (alpha_star_max - alpha_star_min)))
              else
                beta = alpha_star_max / alpha_v
              end if
              if (beta < 0 .or. beta > 1) then
                call log_error('Vertical velocity split weight is out range from 0 to 1!',__FILE__, __LINE__)
              end if
              this%we_imp(i,j,k) = (1 - beta) * this%we(i,j,k)
              this%we    (i,j,k) = beta * this%we(i,j,k)
            end do
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_jds
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            if (this%we(mesh%full_ids,j,k) >= 0) then
              do i = mesh%full_ids, mesh%full_ide
                work(i,k) = max(this%v(i,j,k-1), 0.)
              end do
            else
              do i = mesh%full_ids, mesh%full_ide
                work(i,k) = max(this%v(i,j,k), 0.)
              end do
            end if
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            alpha_h = dt_ * pole(k) * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
            alpha_star_max = alpha_max - kesi * alpha_h
            alpha_star_min = alpha_min * alpha_star_max / alpha_max
            alpha_v = dt_ * abs(this%we(mesh%full_ids,j,k) / this%m(mesh%full_ids,j,k))
            if (alpha_v <= alpha_star_min) then
              beta = 1
            else if (alpha_v > alpha_star_min .and. alpha_v <= 2 * alpha_star_max - alpha_star_min) then
              beta = 1 / (1 + (alpha_v - alpha_star_min)**2 / (4 * alpha_star_max * (alpha_star_max - alpha_star_min)))
            else
              beta = alpha_star_max / alpha_v
            end if
            if (beta < 0 .or. beta > 1) then
              call log_error('Vertial velocity split weight is out of 0 to 1!', __FILE__, __LINE__)
            end if
            this%we_imp(:,j,k) = (1 - beta) * this%we(mesh%full_ids,j,k)
            this%we    (:,j,k) = beta * this%we(mesh%full_ids,j,k)
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_jde
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            if (this%we(mesh%full_ids,j,k) >= 0) then
              do i = mesh%full_ids, mesh%full_ide
                work(i,k) = -min(this%v(i,j,k-1), 0.)
              end do
            else
              do i = mesh%full_ids, mesh%full_ide
                work(i,k) = -min(this%v(i,j,k), 0.)
              end do
            end if
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            alpha_h = dt_ * pole(k) * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
            alpha_star_max = alpha_max - kesi * alpha_h
            alpha_star_min = alpha_min * alpha_star_max / alpha_max
            alpha_v = dt_ * abs(this%we(mesh%full_ids,j,k) / this%m(mesh%full_ids,j,k))
            if (alpha_v <= alpha_star_min) then
              beta = 1
            else if (alpha_v > alpha_star_min .and. alpha_v <= 2 * alpha_star_max - alpha_star_min) then
              beta = 1 / (1 + (alpha_v - alpha_star_min)**2 / (4 * alpha_star_max * (alpha_star_max - alpha_star_min)))
            else
              beta = alpha_star_max / alpha_v
            end if
            if (beta < 0 .or. beta > 1) then
              call log_error('Vertical velocity split weight is out of 0 to 1!', __FILE__, __LINE__)
            end if
            this%we_imp(:,j,k) = (1 - beta) * this%we(mesh%full_ids,j,k)
            this%we    (:,j,k) = beta * this%we(mesh%full_ids,j,k)
          end do
        end if
        this%we    (:,:,mesh%half_kds) = 0; this%we    (:,:,mesh%half_kde) = 0
        this%we_imp(:,:,mesh%half_kds) = 0; this%we_imp(:,:,mesh%half_kde) = 0
      end if
    end if
    end associate

  end subroutine adv_batch_accum_we_lev

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
