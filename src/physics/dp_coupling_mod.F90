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

    integer i, j, k, icol

    associate (mesh     => block%mesh                , &
               pstate   => block%pstate              , & ! out
               u        => block%dstate(itime)%u     , & ! in
               v        => block%dstate(itime)%v     , & ! in
               pt       => block%dstate(itime)%pt    , & ! in
               t        => block%dstate(itime)%t     , & ! in
               p        => block%dstate(itime)%p     , & ! in
               dp       => block%dstate(itime)%m     , & ! in
               p_lev    => block%dstate(itime)%p_lev , & ! in
               gz       => block%dstate(itime)%gz    , & ! in
               gz_lev   => block%dstate(itime)%gz_lev, & ! in
               ps       => block%dstate(itime)%phs   , & ! in
               qv       => block%dstate(itime)%qv    )   ! in
    ! Full levels
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      icol = 0
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          icol = icol + 1
          pstate%u    (icol,k) = u     (i,j,k)
          pstate%v    (icol,k) = v     (i,j,k)
          pstate%pt   (icol,k) = pt    (i,j,k)
          pstate%t    (icol,k) = temperature(pt(i,j,k), p(i,j,k), qv(i,j,k))
          pstate%sh   (icol,k) = specific_humidity(qv(i,j,k))
          pstate%qv   (icol,k) = qv    (i,j,k)
          pstate%p    (icol,k) = p     (i,j,k)
          pstate%p_exn(icol,k) = (p(i,j,k) / p0)**Rd_o_cpd
          pstate%dp   (icol,k) = dp    (i,j,k)
          pstate%z    (icol,k) = gz    (i,j,k) / g
          pstate%rdp  (icol,k) = 1.0_r8 / dp(i,j,k)
          pstate%dz   (icol,k) = (gz_lev(i,j,k+1) - gz_lev(i,j,k)) / g
          pstate%rho  (icol,k) = moist_air_density(t(i,j,k), p(i,j,k), qv(i,j,k))
        end do
      end do
    end do
    ! Half levels
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      icol = 0
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          icol = icol + 1
          pstate%p_lev(icol,k) =  p_lev(i,j,k)
          pstate%z_lev(icol,k) = gz_lev(i,j,k) / g
        end do
      end do
    end do 
    ! Surface
    icol = 0
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
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
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        icol = 0
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            icol = icol + 1
            dudt(i,j,k) = ptend%dudt(icol,k)
            dvdt(i,j,k) = ptend%dvdt(icol,k)
          end do
        end do
      end do
      call fill_halo(block, dudt, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block, dvdt, full_lon=.true., full_lat=.true., full_lev=.true., &
                     west_halo=.false.,  east_halo=.false., south_halo=.false.)
    end if
    if (ptend%updated_t) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        icol = 0
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            icol = icol + 1
            dtdt(i,j,k) = ptend%dtdt(icol,k)
          end do
        end do
      end do
    end if
    if (ptend%updated_sh) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        icol = 0
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            icol = icol + 1
            dshdt(i,j,k) = ptend%dshdt(icol,k)
          end do
        end do
      end do
    end if
    end associate

  end subroutine dp_coupling_p2d

end module dp_coupling_mod
