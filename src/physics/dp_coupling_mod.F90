module dp_coupling_mod

  use const_mod
  use block_mod
  use physics_types_mod
  use formula_mod
  use parallel_mod

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

    associate (mesh     => block%mesh                   , &
               pstate   => block%pstate                 , & ! out
               u        => block%dstate(itime)%u        , & ! in
               v        => block%dstate(itime)%v        , & ! in
               pt       => block%dstate(itime)%pt       , & ! in
               t        => block%dstate(itime)%t        , & ! in
               ph       => block%dstate(itime)%ph       , & ! in
               ph_lev   => block%dstate(itime)%ph_lev   , & ! in
               dph      => block%dstate(itime)%dmg      , & ! in
               p        => block%dstate(itime)%p        , & ! in
               p_lev    => block%dstate(itime)%p_lev    , & ! in
               gz       => block%dstate(itime)%gz       , & ! in
               gz_lev   => block%dstate(itime)%gz_lev   , & ! in
               ps       => block%dstate(itime)%phs      , & ! in
               qv       => block%dstate(itime)%qv       )   ! in
    ! Full levels
    do k = mesh%full_kds, mesh%full_kde
      icol = 0
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          icol = icol + 1
          pstate%u        (icol,k) = u(i,j,k)
          pstate%v        (icol,k) = v(i,j,k)
          pstate%pt       (icol,k) = pt(i,j,k)
          pstate%t        (icol,k) = temperature(pt(i,j,k), p(i,j,k), qv(i,j,k))
          pstate%sh       (icol,k) = specific_humidity(qv(i,j,k))
          pstate%qv       (icol,k) = qv(i,j,k)
          pstate%tv       (icol,k) = virtual_temperature(pstate%t(icol,k), qv(i,j,k))
          pstate%ptv      (icol,k) = virtual_potential_temperature(pstate%tv(icol,k), p(i,j,k))
          pstate%ph       (icol,k) = ph(i,j,k)
          pstate%pkh      (icol,k) = ph(i,j,k)**rd_o_cpd / pk0
          pstate%dph      (icol,k) = dph(i,j,k)
          pstate%p        (icol,k) = p(i,j,k)
          pstate%dp       (icol,k) = p(i,j,k+1) - p(i,j,k)
          pstate%rdp      (icol,k) = 1.0_r8 / pstate%dp(icol,k)
          pstate%z        (icol,k) = gz    (i,j,k) / g
          pstate%dz       (icol,k) = (gz_lev(i,j,k+1) - gz_lev(i,j,k)) / g
          pstate%rho      (icol,k) = moist_air_density(t(i,j,k), p(i,j,k), qv(i,j,k))
          pstate%cp       (icol,k) = (1 - pstate%sh(icol,k)) * cpd + pstate%sh(icol,k) * cpv
          pstate%cv       (icol,k) = (1 - pstate%sh(icol,k)) * cvd + pstate%sh(icol,k) * cvv
          tmp = gz(i,j,k) + 0.5_r8 * (u(i,j,k)**2 + v(i,j,k)**2)
          pstate%tep      (icol,k) = pstate%cp(icol,k) * pstate%t(icol,k) + tmp
          pstate%tev      (icol,k) = pstate%cv(icol,k) * pstate%t(icol,k) + tmp
        end do
      end do
    end do
    ! Half levels
    do k = mesh%half_kds, mesh%half_kde
      icol = 0
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          icol = icol + 1
          pstate%ph_lev   (icol,k) = ph_lev(i,j,k)
          pstate%lnph_lev (icol,k) = log(ph_lev(i,j,k))
          pstate%z_lev    (icol,k) = gz_lev(i,j,k) / g
          if (mesh%half_kds < k .and. k < mesh%half_kde) then
            pstate%n2_lev(icol,k) = buoyancy_frequency( &
              pstate%ptv(icol,k-1), pstate%ptv(icol,k), pstate%z(icol,k-1), pstate%z(icol,k))
            pstate%ri_lev(icol,k) = local_richardson_number( &
              pstate%N2_lev(icol,k), pstate%z(icol,k-1), pstate%z(icol,k), u(i,j,k-1), u(i,j,k), v(i,j,k-1), v(i,j,k))
          end if
        end do
      end do
    end do 
    ! Surface
    icol = 0
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        icol = icol + 1
        pstate%ps(icol) = ps(i,j)
      end do
    end do
    end associate

  end subroutine dp_coupling_d2p

  subroutine dp_coupling_p2d(block, itime)

    type(block_type), intent(inout) :: block    
    integer, intent(in) :: itime

    integer i, j, k, icol

    associate (mesh  => block%mesh                   , &
               ptend => block%ptend                  , & ! in
               dudt  => block%dtend(itime)%dudt_phys , & ! out
               dvdt  => block%dtend(itime)%dvdt_phys , & ! out
               dtdt  => block%dtend(itime)%dtdt_phys , & ! out
               dshdt => block%dtend(itime)%dshdt_phys)   ! out
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
      call fill_halo(block%halo, dudt, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block%halo, dvdt, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false.,  east_halo=.false., south_halo=.false.)
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
    if (ptend%updated_sh) then
      do k = mesh%full_kds, mesh%full_kde
        icol = 0
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            icol = icol + 1
            dshdt(i,j,k) = ptend%dshdt(icol,k)
          end do
        end do
      end do
    end if
    end associate

  end subroutine dp_coupling_p2d

end module dp_coupling_mod
