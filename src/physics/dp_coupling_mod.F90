! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module dp_coupling_mod

  use const_mod
  use namelist_mod
  use block_mod
  use physics_types_mod
  use formula_mod
  use latlon_parallel_mod
  use tracer_mod
#ifdef HAS_CAM
  use cam_physics_driver_mod, only: cam_physics_d2p, cam_physics_p2d
#endif
  use filter_mod

  implicit none

  private

  public dp_coupling_d2p
  public dp_coupling_p2d

contains

  subroutine dp_coupling_d2p(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    real(r8) tmp
    integer i, j, k, icol

    if (physics_suite == 'N/A') return

    select case (physics_suite)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_d2p(block, itime)
#endif
    case default
      associate (mesh   => block%mesh                , &
                 pstate => block%pstate              , & ! out
                 u      => block%dstate(itime)%u     , & ! in
                 v      => block%dstate(itime)%v     , & ! in
                 pt     => block%dstate(itime)%pt    , & ! in (modified potential temperature)
                 t      => block%dstate(itime)%t     , & ! in
                 tv     => block%dstate(itime)%tv    , & ! in (virtual temperature)
                 q      => tracers(block%id)%q       , & ! in (tracer dry mixing ratio)
                 qm     => tracers(block%id)%qm      , & ! in (total moist dry mixing ratio)
                 ph     => block%dstate(itime)%ph    , & ! in (hydrostatic full pressure)
                 ph_lev => block%dstate(itime)%ph_lev, & ! in
                 phs    => block%dstate(itime)%phs   , & ! in (surface hydrostatic pressure)
                 dmg    => block%dstate(itime)%dmg   , & ! in (dry air weight within each layer)
                 omg    => block%aux%omg             , & ! in (vertical pressure velocity)
                 gz     => block%dstate(itime)%gz    , & ! in (geopotential)
                 gz_lev => block%dstate(itime)%gz_lev, & ! in
                 land   => block%static%landmask     )   ! in
      ! Full levels
      do k = mesh%full_kds, mesh%full_kde
        icol = 0
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            icol = icol + 1
            pstate%u        (icol,k)   = u %d(i,j,k)
            pstate%v        (icol,k)   = v %d(i,j,k)
            pstate%pt       (icol,k)   = pt%d(i,j,k) ! FIXME: What does physics need? Dry air potential temperature or just modified one?
            pstate%t        (icol,k)   = t %d(i,j,k)
            pstate%tv       (icol,k)   = tv%d(i,j,k)
            pstate%ptv      (icol,k)   = virtual_potential_temperature(tv%d(i,j,k), ph%d(i,j,k))
            pstate%q        (icol,k,:) = wet_mixing_ratio(q%d(i,j,k,:), qm%d(i,j,k))
            pstate%p        (icol,k)   = ph    %d(i,j,k)
            pstate%p_lev    (icol,k)   = ph_lev%d(i,j,k)
            pstate%pk       (icol,k)   = ph    %d(i,j,k)**rd_o_cpd / pk0
            pstate%pk_lev   (icol,k)   = ph_lev%d(i,j,k)**rd_o_cpd / pk0
            pstate%dp       (icol,k)   = ph_lev%d(i,j,k+1) - ph_lev%d(i,j,k)
            pstate%rdp      (icol,k)   = 1.0_r8 / pstate%dp(icol,k)
            pstate%omg      (icol,k)   = omg%d(i,j,k)
            pstate%z        (icol,k)   = gz%d(i,j,k) / g
            pstate%dz       (icol,k)   = (gz_lev%d(i,j,k+1) - gz_lev%d(i,j,k)) / g
            pstate%rho      (icol,k)   = moist_air_density(t%d(i,j,k), ph%d(i,j,k), pstate%qv(icol,k), qm%d(i,j,k))
            pstate%cp       (icol,k)   = (1 - pstate%qv(icol,k)) * cpd + pstate%qv(icol,k) * cpv ! FIXME: Add specific heat capacities of other water species.
            pstate%cv       (icol,k)   = (1 - pstate%qv(icol,k)) * cvd + pstate%qv(icol,k) * cvv ! FIXME: Add specific heat capacities of other water species.
            tmp = gz%d(i,j,k) + 0.5_r8 * (u%d(i,j,k)**2 + v%d(i,j,k)**2)
            pstate%tep      (icol,k)   = pstate%cp(icol,k) * pstate%t(icol,k) + tmp
            pstate%tev      (icol,k)   = pstate%cv(icol,k) * pstate%t(icol,k) + tmp
          end do
        end do
      end do
      ! Half levels
      do k = mesh%half_kds, mesh%half_kde
        icol = 0
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            icol = icol + 1
            pstate%p_lev   (icol,k) = ph_lev%d(i,j,k)
            pstate%lnp_lev (icol,k) = log(ph_lev%d(i,j,k))
            pstate%z_lev   (icol,k) = gz_lev%d(i,j,k) / g
            if (mesh%half_kds < k .and. k < mesh%half_kde) then
              pstate%n2_lev(icol,k) = buoyancy_frequency( &
                pstate%ptv   (icol,k-1), &
                pstate%ptv   (icol,k  ), &
                pstate%z     (icol,k-1), &
                pstate%z     (icol,k  ))
              pstate%ri_lev(icol,k) = local_richardson_number( &
                pstate%N2_lev(icol,k  ), &
                pstate%z     (icol,k-1), &
                pstate%z     (icol,k  ), &
                u          %d(i,j ,k-1), &
                u          %d(i,j ,k  ), &
                v          %d(i,j ,k-1), &
                v          %d(i,j ,k  ))
            end if
          end do
        end do
      end do
      ! Surface
      icol = 0
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          icol = icol + 1
          pstate%ps  (icol) = phs%d(i,j)
          pstate%wsb (icol) = sqrt(u%d(i,j,mesh%full_kde)**2 + v%d(i,j,mesh%full_kde)**2)
          pstate%land(icol) = land%d(i,j)
        end do
      end do
      end associate
    end select

  end subroutine dp_coupling_d2p

  subroutine dp_coupling_p2d(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer i, j, k, icol, m

    select case (physics_suite)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_p2d(block)
#endif
    case default
      associate (mesh  => block%mesh             , &
                 ptend => block%ptend            , & ! in
                 dmg   => block%dstate(itime)%dmg, & ! in
                 ph    => block%dstate(itime)%ph , & ! in
                 t     => block%dstate(itime)%t  , & ! in
                 q     => tracers(block%id)%q    , & ! in (tracer dry mixing ratio)
                 dudt  => block%aux%dudt_phys    , & ! out
                 dvdt  => block%aux%dvdt_phys    , & ! out
                 dptdt => block%aux%dptdt_phys   , & ! out
                 dqdt  => block%aux%dqdt_phys)       ! out
      if (ptend%updated_u .and. ptend%updated_v) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 0
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              icol = icol + 1
              dudt%d(i,j,k) = ptend%dudt(icol,k)
              dvdt%d(i,j,k) = ptend%dvdt(icol,k)
            end do
          end do
        end do
      end if
      do m = 1, ntracers
        if (ptend%updated_q(m)) then
          do k = mesh%full_kds, mesh%full_kde
            icol = 0
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                icol = icol + 1
                ! Convert to dry mixing ratio tendency.
                dqdt%d(i,j,k,m) = ptend%dqdt(icol,k,m) / (1 - q%d(i,j,k,m))**2
              end do
            end do
          end do
        end if
      end do
      if (ptend%updated_t) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 0
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              icol = icol + 1
              ! Convert to modified potential temperature tendency and multiply dmg.
              dptdt%d(i,j,k) = dmg%d(i,j,k) * (p0 / ph%d(i,j,k))**rd_o_cpd * ( &
                (1 + rv_o_rd * q%d(i,j,k,idx_qv)) * ptend%dtdt(icol,k) + &
                rv_o_rd * t%d(i,j,k) * dqdt%d(i,j,k,idx_qv))
            end do
          end do
        end do
      end if
      do m = 1, ntracers
        if (ptend%updated_q(m)) then
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                ! Multiply dmg.
                dqdt%d(i,j,k,m) = dmg%d(i,j,k) * dqdt%d(i,j,k,m)
              end do
            end do
          end do
        end if
      end do
      end associate
    end select

    if (filter_ptend) then
      call fill_halo(block%aux%dudt_phys, south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dudt_phys%d)
      call fill_halo(block%aux%dvdt_phys, south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dvdt_phys%d)
      call fill_halo(block%aux%dptdt_phys, south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dptdt_phys%d)
      do m = 1, ntracers
        call fill_halo(block%aux%dqdt_phys, m, south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, block%aux%dqdt_phys%d(:,:,:,m))
      end do
      call fill_halo(block%aux%dudt_phys, west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block%aux%dvdt_phys, west_halo=.false.,  east_halo=.false., south_halo=.false.)
    else
      call fill_halo(block%aux%dudt_phys, west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block%aux%dvdt_phys, west_halo=.false.,  east_halo=.false., south_halo=.false.)
    end if

  end subroutine dp_coupling_p2d

end module dp_coupling_mod
