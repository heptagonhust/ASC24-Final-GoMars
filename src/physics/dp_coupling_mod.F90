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

    real(r8), pointer :: q(:,:,:,:), qm(:,:,:)
    real(r8) tmp
    integer i, j, k, icol

    if (physics_suite == 'N/A') return

    select case (physics_suite)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_d2p(block, itime)
#endif
    case default
      call tracer_get_array(block%id, q)
      call tracer_get_array_qm(block%id, qm)
      associate (mesh        => block%mesh                   , &
                 pstate      => block%pstate                 , & ! out
                 u           => block%dstate(itime)%u        , & ! in
                 v           => block%dstate(itime)%v        , & ! in
                 pt          => block%dstate(itime)%pt       , & ! in (modified potential temperature)
                 t           => block%dstate(itime)%t        , & ! in
                 tv          => block%dstate(itime)%tv       , & ! in (virtual temperature)
                 p           => block%dstate(itime)%ph       , & ! in (hydrostatic full pressure)
                 p_lev       => block%dstate(itime)%ph_lev   , & ! in
                 ps          => block%dstate(itime)%phs      , & ! in (surface hydrostatic pressure)
                 dmg         => block%dstate(itime)%dmg      , & ! in (dry air weight within each layer)
                 omg         => block%aux%omg                , & ! in (vertical pressure velocity)
                 gz          => block%dstate(itime)%gz       , & ! in (geopotential)
                 gz_lev      => block%dstate(itime)%gz_lev   , & ! in
                 land        => block%static%landmask        )   ! in
      ! Full levels
      do k = mesh%full_kds, mesh%full_kde
        icol = 0
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            icol = icol + 1
            pstate%u        (icol,k)   = u (i,j,k)
            pstate%v        (icol,k)   = v (i,j,k)
            pstate%pt       (icol,k)   = pt(i,j,k) ! FIXME: What does physics need? Dry air potential temperature or just modified one?
            pstate%t        (icol,k)   = t (i,j,k)
            pstate%tv       (icol,k)   = tv(i,j,k)
            pstate%ptv      (icol,k)   = virtual_potential_temperature(tv(i,j,k), p(i,j,k))
            pstate%q        (icol,k,:) = wet_mixing_ratio(q(i,j,k,:), qm(i,j,k))
            pstate%p        (icol,k)   = p    (i,j,k)
            pstate%p_lev    (icol,k)   = p_lev(i,j,k)
            pstate%pk       (icol,k)   = p    (i,j,k)**rd_o_cpd / pk0
            pstate%pk_lev   (icol,k)   = p_lev(i,j,k)**rd_o_cpd / pk0
            pstate%dp       (icol,k)   = p_lev(i,j,k+1) - p_lev(i,j,k)
            pstate%rdp      (icol,k)   = 1.0_r8 / pstate%dp(icol,k)
            pstate%omg      (icol,k)   = omg(i,j,k)
            pstate%z        (icol,k)   = gz(i,j,k) / g
            pstate%dz       (icol,k)   = (gz_lev(i,j,k+1) - gz_lev(i,j,k)) / g
            pstate%rho      (icol,k)   = moist_air_density(t(i,j,k), p(i,j,k), pstate%qv(icol,k), qm(i,j,k))
            pstate%cp       (icol,k)   = (1 - pstate%qv(icol,k)) * cpd + pstate%qv(icol,k) * cpv ! FIXME: Add specific heat capacities of other water species.
            pstate%cv       (icol,k)   = (1 - pstate%qv(icol,k)) * cvd + pstate%qv(icol,k) * cvv ! FIXME: Add specific heat capacities of other water species.
            tmp = gz(i,j,k) + 0.5_r8 * (u(i,j,k)**2 + v(i,j,k)**2)
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
            pstate%p_lev   (icol,k) = p_lev(i,j,k)
            pstate%lnp_lev (icol,k) = log(p_lev(i,j,k))
            pstate%z_lev   (icol,k) = gz_lev(i,j,k) / g
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
                u            (i,j ,k-1), &
                u            (i,j ,k  ), &
                v            (i,j ,k-1), &
                v            (i,j ,k  ))
            end if
          end do
        end do
      end do
      ! Surface
      icol = 0
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          icol = icol + 1
          pstate%ps  (icol) = ps(i,j)
          pstate%wsb (icol) = sqrt(u(i,j,mesh%full_kde)**2 + v(i,j,mesh%full_kde)**2)
          pstate%land(icol) = land(i,j)
        end do
      end do
      end associate
    end select

  end subroutine dp_coupling_d2p

  subroutine dp_coupling_p2d(block)

    type(block_type), intent(inout) :: block

    integer i, j, k, icol, m

    select case (physics_suite)
#ifdef HAS_CAM
    case ('cam')
      call cam_physics_p2d(block)
#endif
    case default
      associate (mesh  => block%mesh          , &
                 ptend => block%ptend         , & ! in
                 dudt  => block%aux%dudt_phys , & ! out
                 dvdt  => block%aux%dvdt_phys , & ! out
                 dtdt  => block%aux%dtdt_phys , & ! out
                 dqdt  => block%aux%dqdt_phys)    ! out
      if (ptend%updated_u .and. ptend%updated_v) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 0
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              icol = icol + 1
              dudt(i,j,k) = ptend%dudt(icol,k)
              dvdt(i,j,k) = ptend%dvdt(icol,k)
            end do
          end do
        end do
      end if
      if (ptend%updated_t) then
        do k = mesh%full_kds, mesh%full_kde
          icol = 0
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              icol = icol + 1
              dtdt(i,j,k) = ptend%dtdt(icol,k)
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
                dqdt(i,j,k,m) = ptend%dqdt(icol,k,m)
              end do
            end do
          end do
        end if
      end do
      end associate
    end select

    if (filter_ptend) then
      call fill_halo(block%filter_halo, block%aux%dudt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dudt_phys)
      call fill_halo(block%filter_halo, block%aux%dvdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dvdt_phys)
      call fill_halo(block%filter_halo, block%aux%dtdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%big_filter, block%aux%dtdt_phys)
      do m = 1, ntracers
        call fill_halo(block%filter_halo, block%aux%dqdt_phys(:,:,:,m), full_lon=.true., full_lat=.true., full_lev=.true., &
                       south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, block%aux%dqdt_phys(:,:,:,m))
      end do
      call fill_halo(block%filter_halo, block%aux%dudt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block%filter_halo, block%aux%dvdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false.,  east_halo=.false., south_halo=.false.)
    else
      call fill_halo(block%halo, block%aux%dudt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block%halo, block%aux%dvdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false.,  east_halo=.false., south_halo=.false.)
    end if

  end subroutine dp_coupling_p2d

end module dp_coupling_mod
