module physics_mod

  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use tracer_mod
  use physics_types_mod
  use dp_coupling_mod
  use latlon_parallel_mod
  use pbl_driver_mod
#ifdef HAS_CCPP
  use ccpp_driver_mod
#endif
#ifdef HAS_CAM
  use cam_physics_driver_mod
#endif

  implicit none

  private

  public physics_init
  public physics_run
  public physics_update
  public physics_update_dynamics
  public physics_update_tracers
  public physics_final

contains

  subroutine physics_init(namelist_path)

    character(*), intent(in) :: namelist_path

    call time_add_alert('phys', seconds=dt_phys)

    select case (physics_suite)
    case ('ccpp')
#ifdef HAS_CCPP
      call ccpp_driver_init(namelist_path)
#else
      if (proc%is_root()) call log_error('CCPP physics is not compiled!')
#endif
    case ('cam')
#ifdef HAS_CAM
      call cam_physics_driver_init(namelist_path)
#else
      if (proc%is_root()) call log_error('CAM physics is not compiled!')
#endif
    end select

  end subroutine physics_init

  subroutine physics_run(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    if (.not. time_is_alerted('phys')) return

    associate (mesh   => block%mesh  , &
               pstate => block%pstate, &
               ptend  => block%ptend )
    call dp_coupling_d2p(block, itime)

    call ptend%reset()

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics(      &
        pstate%ncol           , &
        pstate%nlev           , &
        dt                    , &
        pstate%lat * rad      , &
        pstate%t              , &
        pstate%qv             , &
        pstate%u              , &
        pstate%v              , &
        pstate%p              , &
        pstate%p_lev          , &
        pstate%dp             , &
        pstate%rdp            , &
        pstate%ps             , &
        ptend%dudt            , &
        ptend%dvdt            , &
        ptend%dtdt            , &
        ptend%dqdt(:,:,idx_qv), &
        pstate%precl          , &
        0                     , & ! test
        .true.                , & ! RJ2012_precip
        .false.                 & ! TC_PBL_mod
      )
      ptend%updated_u    = .true.
      ptend%updated_v    = .true.
      ptend%updated_t    = .true.
      ptend%updated_q(1) = .true.
#ifdef HAS_CCPP
    case ('ccpp')
      call ccpp_driver_run()
#endif
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_driver_run1()
      call cam_physics_driver_sfc_flux()
      call cam_physics_driver_run2()
#endif
    end select

    call dp_coupling_p2d(block, itime)
    end associate

  end subroutine physics_run

  subroutine physics_update(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    real(r8), pointer :: q(:,:,:,:)
    integer i, j, k, m

    associate (mesh  => block%mesh               , &
               ptend => block%ptend              , &
               aux   => block%aux                , &
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   )   ! inout
    ! Update dynamics.
    if (ptend%updated_u) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u_lon(i,j,k) = u_lon(i,j,k) + dt * 0.5_r8 * (aux%dudt_phys(i,j,k) + aux%dudt_phys(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true. , full_lev=.true.)
    end if
    if (ptend%updated_v) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v_lat(i,j,k) = v_lat(i,j,k) + dt * 0.5_r8 * (aux%dvdt_phys(i,j,k) + aux%dvdt_phys(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(block%halo, v_lat, full_lon=.true. , full_lat=.false., full_lev=.true.)
    end if
    if (ptend%updated_t) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt(i,j,k) = pt(i,j,k) + dt * aux%dptdt_phys(i,j,k) / dmg(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
    end if
    ! Update tracers.
    call tracer_get_array(block%id, q)
    do m = 1, ntracers
      if (ptend%updated_q(m)) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q(i,j,k,m) = q(i,j,k,m) + dt * aux%dqdt_phys(i,j,k,m) / dmg(i,j,k)
            end do
          end do
        end do
        call tracer_fill_negative_values(block, itime, q(:,:,:,m))
        call fill_halo(block%filter_halo, q(:,:,:,m), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
      end if
    end do
    call tracer_calc_qm(block)
    end associate

  end subroutine physics_update

  subroutine physics_update_dynamics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh  => block%mesh               , &
               ptend => block%ptend              , &
               aux   => block%aux                , &
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   )   ! inout
    ! Update dynamics.
    if (ptend%updated_u) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u_lon(i,j,k) = u_lon(i,j,k) + dt * 0.5_r8 * (aux%dudt_phys(i,j,k) + aux%dudt_phys(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true. , full_lev=.true.)
    end if
    if (ptend%updated_v) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v_lat(i,j,k) = v_lat(i,j,k) + dt * 0.5_r8 * (aux%dvdt_phys(i,j,k) + aux%dvdt_phys(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(block%halo, v_lat, full_lon=.true. , full_lat=.false., full_lev=.true.)
    end if
    if (ptend%updated_t) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt(i,j,k) = pt(i,j,k) + dt * aux%dptdt_phys(i,j,k) / dmg(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
    end if
    end associate

  end subroutine physics_update_dynamics

  subroutine physics_update_tracers(block, itime, dt, idx)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt
    integer, intent(in) :: idx

    real(r8), pointer :: q(:,:,:)
    integer i, j, k, m

    associate (mesh  => block%mesh             , &
               ptend => block%ptend            , &
               aux   => block%aux              , &
               dmg   => block%dstate(itime)%dmg)
    ! Update tracers.
    call tracer_get_array(block%id, idx, q, __FILE__, __LINE__)
    if (ptend%updated_q(idx)) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            q(i,j,k) = q(i,j,k) + dt * aux%dqdt_phys(i,j,k,idx) / dmg(i,j,k)
          end do
        end do
      end do
      call tracer_fill_negative_values(block, itime, q)
      call fill_halo(block%filter_halo, q, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
    end if
    end associate

  end subroutine physics_update_tracers

  subroutine physics_final()

    select case (physics_suite)
#ifdef HAS_CCPP
    case ('ccpp')
      call ccpp_driver_final()
#endif
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_driver_final()
#endif
    end select

  end subroutine physics_final

end module physics_mod
