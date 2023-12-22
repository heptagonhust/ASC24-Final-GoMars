! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

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
  use mars_nasa_physics_driver_mod

  implicit none

  private

  public physics_init_stage1
  public physics_init_stage2
  public physics_run
  public physics_update
  public physics_update_dynamics
  public physics_update_tracers
  public physics_final

contains

  subroutine physics_init_stage1(namelist_path)

    character(*), intent(in) :: namelist_path

    if (physics_suite /= 'N/A') then
      call time_add_alert('phys', seconds=dt_phys)
    end if

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
    case ('mars_nasa')
      call mars_nasa_physics_driver_init(namelist_path, nlev)
    end select

  end subroutine physics_init_stage1

  subroutine physics_init_stage2()


  end subroutine physics_init_stage2

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
      ptend%updated_u         = .true.
      ptend%updated_v         = .true.
      ptend%updated_t         = .true.
      ptend%updated_q(idx_qv) = .true.
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

    integer i, j, k, m

    associate (mesh  => block%mesh               , &
               ptend => block%ptend              , &
               dudt  => block%aux%dudt_phys      , & ! in
               dvdt  => block%aux%dvdt_phys      , & ! in
               dptdt => block%aux%dptdt_phys     , & ! in
               dqdt  => block%aux%dqdt_phys      , & ! in
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   , & ! inout
               q     => tracers(block%id)%q      )   ! inout
    ! Update dynamics.
    if (ptend%updated_u) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * 0.5_r8 * (dudt%d(i,j,k) + dudt%d(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(u_lon)
    end if
    if (ptend%updated_v) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * 0.5_r8 * (dvdt%d(i,j,k) + dvdt%d(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(v_lat)
    end if
    if (ptend%updated_t) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt%d(i,j,k) / dmg%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(pt)
    end if
    ! Update tracers.
    do m = 1, ntracers
      if (ptend%updated_q(m)) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,m) = q%d(i,j,k,m) + dt * dqdt%d(i,j,k,m) / dmg%d(i,j,k)
            end do
          end do
        end do
        call tracer_fill_negative_values(block, itime, q%d(:,:,:,m))
        call fill_halo(q, m)
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
               dudt  => block%aux%dudt_phys      , & ! in
               dvdt  => block%aux%dvdt_phys      , & ! in
               dptdt => block%aux%dptdt_phys     , & ! in
               dmg   => block%dstate(itime)%dmg  , & ! in
               u_lon => block%dstate(itime)%u_lon, & ! inout
               v_lat => block%dstate(itime)%v_lat, & ! inout
               pt    => block%dstate(itime)%pt   )   ! inout
    ! Update dynamics.
    if (ptend%updated_u) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u_lon%d(i,j,k) = u_lon%d(i,j,k) + dt * 0.5_r8 * (dudt%d(i,j,k) + dudt%d(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(u_lon)
    end if
    if (ptend%updated_v) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v_lat%d(i,j,k) = v_lat%d(i,j,k) + dt * 0.5_r8 * (dvdt%d(i,j,k) + dvdt%d(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(v_lat)
    end if
    if (ptend%updated_t) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt%d(i,j,k) = pt%d(i,j,k) + dt * dptdt%d(i,j,k) / dmg%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(pt)
    end if
    end associate

  end subroutine physics_update_dynamics

  subroutine physics_update_tracers(block, itime, dt, idx)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt
    integer, intent(in) :: idx

    integer i, j, k

    associate (mesh  => block%mesh             , &
               ptend => block%ptend            , &
               dqdt  => block%aux%dqdt_phys    , &
               dmg   => block%dstate(itime)%dmg, &
               q     => tracers(block%id)%q    )
    ! Update tracers.
    if (ptend%updated_q(idx)) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            q%d(i,j,k,idx) = q%d(i,j,k,idx) + dt * dqdt%d(i,j,k,idx) / dmg%d(i,j,k)
          end do
        end do
      end do
      call tracer_fill_negative_values(block, itime, q%d(:,:,:,idx))
      call fill_halo(q, idx)
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
