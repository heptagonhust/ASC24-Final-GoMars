module physics_mod

  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use physics_types_mod
  use moist_mod
  use dp_coupling_mod
  use parallel_mod
  use formula_mod
  use pbl_driver_mod

  implicit none

  private

  public physics_init
  public physics_run_before_dynamics
  public physics_run_after_dynamics
  public physics_final

contains

  subroutine physics_init()

    call time_add_alert('phys', seconds=dt_phys)

  end subroutine physics_init

  subroutine physics_run_before_dynamics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(8), intent(in) :: dt

  end subroutine physics_run_before_dynamics

  subroutine physics_run_after_dynamics(block, itime, dt)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime
    real(8), intent(in) :: dt

    integer i, j, k

    if (.not. time_is_alerted('phys')) return

    associate (mesh   => block%mesh         , &
               pstate => block%pstate       , &
               ptend  => block%ptend        , &
               dstate => block%dstate(itime), &
               dtend   => block%dtend(itime))
    call dp_coupling_d2p(block, itime)

    call ptend%reset()

    select case (physics_suite)
    case ('simple_physics')
      call simple_physics( &
        pstate%ncol , &
        pstate%nlev , &
        dt          , &
        pstate%lat  , &
        pstate%t    , &
        pstate%sh   , &
        pstate%u    , &
        pstate%v    , &
        pstate%p    , &
        pstate%p_lev, &
        pstate%dp   , &
        pstate%rdp  , &
        pstate%ps   , &
        ptend%dudt  , &
        ptend%dvdt  , &
        ptend%dtdt  , &
        ptend%dshdt , &
        pstate%precl, &
        0           , & ! test
        .true.      , & ! RJ2012_precip
        .false.       & ! TC_PBL_mod
      )
      ptend%updated_u  = .true.
      ptend%updated_v  = .true.
      ptend%updated_t  = .true.
      ptend%updated_sh = .true.
    end select

    call dp_coupling_p2d(block, itime)

    if (ptend%updated_u .and. ptend%updated_v) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            dstate%u_lon(i,j,k) = dstate%u_lon(i,j,k) + dt * 0.5_r8 * (dtend%dudt_phys(i,j,k) + dtend%dudt_phys(i+1,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%u_lon, full_lon=.false., full_lat=.true. , full_lev=.true.)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dstate%v_lat(i,j,k) = dstate%v_lat(i,j,k) + dt * 0.5_r8 * (dtend%dvdt_phys(i,j,k) + dtend%dvdt_phys(i,j+1,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%v_lat, full_lon=.true. , full_lat=.false., full_lev=.true.)
    end if

    if (ptend%updated_sh) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dstate%qv(i,j,k) = mixing_ratio(specific_humidity(dstate%qv(i,j,k)) + dt * dtend%dshdt_phys(i,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%qv, full_lon=.true. , full_lat=.true. , full_lev=.true.)
      call calc_qm(block, itime)
    end if

    if (ptend%updated_t) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dstate%t (i,j,k) = dstate%t(i,j,k) + dt * dtend%dtdt_phys(i,j,k)
            dstate%pt(i,j,k) = potential_temperature(dstate%t(i,j,k), dstate%p(i,j,k), dstate%qv(i,j,k))
          end do
        end do
      end do
      call fill_halo(block%halo, dstate%pt, full_lon=.true. , full_lat=.true. , full_lev=.true.)
    end if
    end associate

  end subroutine physics_run_after_dynamics

  subroutine physics_final()

  end subroutine physics_final

end module physics_mod
