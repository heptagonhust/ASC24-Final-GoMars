module steady_state_pgf_test_mod

  use flogger
  use const_mod, only: r8, pi, Rd, g, omega
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use formula_mod
  use operators_mod

  implicit none

  private

  public steady_state_pgf_test_set_params
  public steady_state_pgf_test_set_ic

  real(r8), parameter :: T0    = 300.0d0     ! K
  real(r8), parameter :: h0    = 2000.0d0    ! m
  real(r8), parameter :: p0    = 1.0d5       ! Pa
  real(r8), parameter :: ztop  = 12.0e3      ! m
  real(r8), parameter :: lonc  = 3 * pi / 2
  real(r8), parameter :: latc  = 0
  real(r8), parameter :: Rm    = 3 * pi / 4
  real(r8), parameter :: gamma = 0.0065d0
  real(r8), parameter :: osm   = pi / 16.0d0

contains

  subroutine steady_state_pgf_test_set_params()

    omega = 0.0
    ptop  = p0 * (1 - gamma / T0 * ztop)**(g / rd / gamma)

  end subroutine steady_state_pgf_test_set_params

  subroutine steady_state_pgf_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r, height
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               mg     => block%dstate(1)%mg    , &
               t      => block%dstate(1)%t     , &
               pt     => block%dstate(1)%pt    , &
               gzs    => block%static%gzs)
    u%d = 0
    v%d = 0

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        full_lon = mesh%full_lon(i)
        r = acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        if (r < Rm) gzs%d(i,j) = g * h0 / 2 * (1 + cos(pi * r / Rm)) * cos(pi * r / osm)**2
      end do
    end do
    call fill_halo(gzs)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        mgs%d(i,j) = p0 * (1 - gamma / T0 * gzs%d(i,j) / g)**(g / Rd / gamma)
      end do
    end do
    call fill_halo(mgs)

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          t %d(i,j,k) = T0 * (mg%d(i,j,k) / p0)**(Rd * gamma / g)
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine steady_state_pgf_test_set_ic

end module steady_state_pgf_test_mod
