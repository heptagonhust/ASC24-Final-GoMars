module physics_mod

  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use tracer_mod
  use physics_types_mod
  use dp_coupling_mod
  use latlon_parallel_mod
  use formula_mod
  use operators_mod
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
  public physics_update_state
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

    call dp_coupling_p2d(block)
    end associate

  end subroutine physics_run

  subroutine physics_update_state(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    real(r8), pointer :: q(:,:,:,:), qm(:,:,:)
    integer i, j, k, m

    associate (mesh   => block%mesh , &
               ptend  => block%ptend, &
               aux    => block%aux  )
    if (ptend%updated_u) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            dstate%u_lon(i,j,k) = dstate%u_lon(i,j,k) + dt * 0.5_r8 * (aux%dudt_phys(i,j,k) + aux%dudt_phys(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%u_lon, full_lon=.false., full_lat=.true. , full_lev=.true.)
    end if
    if (ptend%updated_v) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dstate%v_lat(i,j,k) = dstate%v_lat(i,j,k) + dt * 0.5_r8 * (aux%dvdt_phys(i,j,k) + aux%dvdt_phys(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%v_lat, full_lon=.true. , full_lat=.false., full_lev=.true.)
    end if

    ! Update tracers.
    call tracer_get_array(block%id, q)
    call tracer_get_array_qm(block%id, qm)
    do m = 1, ntracers
      if (ptend%updated_q(m)) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q(i,j,k,m) = wet_mixing_ratio(q(i,j,k,m), qm(i,j,k)) + dt * aux%dqdt_phys(i,j,k,m)
            end do
          end do
        end do
      end if
    end do
    call tracer_calc_qm(block)
    do m = 1, ntracers
      if (ptend%updated_q(m)) then
        ! FIXME: When q is not updated, we should only update its value with new qm.
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              q(i,j,k,m) = dry_mixing_ratio(q(i,j,k,m), qm(i,j,k))
            end do
          end do
        end do
        call fill_halo(block%filter_halo, q(:,:,:,m), full_lon=.true. , full_lat=.true. , full_lev=.true., cross_pole=.true.)
      end if
    end do

    if (ptend%updated_t) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dstate%t (i,j,k) = dstate%t(i,j,k) + dt * aux%dtdt_phys(i,j,k)
            dstate%pt(i,j,k) = modified_potential_temperature(dstate%t(i,j,k), dstate%p(i,j,k), q(i,j,k,idx_qv))
          end do
        end do
      end do
      call fill_halo(block%filter_halo, dstate%pt, full_lon=.true. , full_lat=.true. , full_lev=.true., cross_pole=.true.)
    end if
    end associate

  end subroutine physics_update_state

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
