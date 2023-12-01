module rossby_haurwitz_wave_3d_test_mod

  use flogger
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use formula_mod
  use operators_mod

  implicit none

  private

  public rossby_haurwitz_wave_3d_test_set_ic

  integer , parameter :: n     = 4
  real(r8), parameter :: u0    = 50d0               ! m s-1
  real(r8)            :: gz0
  real(r8)            :: M
  real(r8), parameter :: t0    = 288d0              ! K
  real(r8), parameter :: gamma = 0.0065d0           ! K m-1
  real(r8), parameter :: pref  = 955d2              ! Pa

contains

  subroutine rossby_haurwitz_wave_3d_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) lon, cos_lat, sin_lat
    real(r8) a, a_lat, b_lat, c_lat, phi_p
    integer i, j, k

    gz0 = 8.0d3 * g
    M   = u0 / (n * radius)

    block%static%gzs%d = 0d0

    a = radius

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               mgs    => block%dstate(1)%mgs   , &
               pt     => block%dstate(1)%pt    , &
               t      => block%dstate(1)%t     , &
               mg_lev => block%dstate(1)%mg_lev, &
               mg     => block%dstate(1)%mg    , &
               gz_lev => block%dstate(1)%gz_lev, &
               gz     => block%dstate(1)%gz)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        cos_lat = mesh%full_cos_lat(j)
        sin_lat = mesh%full_sin_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          lon = mesh%half_lon(i)
          u%d(i,j,k) = a * M * cos_lat + a * M * cos_lat**(n - 1) * cos(n * lon) * (n * sin_lat**2 - cos_lat**2)
        end do
      end do
    end do
    call fill_halo(u)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        cos_lat = mesh%half_cos_lat(j)
        sin_lat = mesh%half_sin_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          v%d(i,j,k) = - a * M * n * cos_lat**(n-1) * sin_lat * sin(n * lon)
        end do
      end do
    end do
    call fill_halo(v)

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        lon = mesh%full_lon(i)
        a_lat = M * (2 * omega + M) / 2 * cos_lat**2 + M**2 / 4d0 * cos_lat**(2 * n) * ((n + 1) * cos_lat**2 + (2 * n**2 - n - 2)) - n**2 * M**2 / 2 * cos_lat**(2 * (n - 1))
        b_lat = 2 * (omega + M) * M / (n + 1) / (n + 2) * cos_lat**n * ((n**2 + 2 * n + 2) - (n + 1)**2 * cos_lat**2)
        c_lat = M**2 / 4d0 * cos_lat**(2 * n) * ((n + 1) * cos_lat**2 - (n + 2))
        phi_p = a**2 * (a_lat + b_lat * cos(n * lon) + c_lat * cos(2 * n * lon))
        mgs%d(i,j) = pref * (1 + gamma / g / t0 * phi_p)**(g / gamma / Rd)
      end do
    end do
    call fill_halo(mgs)

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          t %d(i,j,k) = t0 * (mg%d(i,j,k) / pref)**(gamma * Rd / g)
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(pt, cross_pole=.true.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gz%d(i,j,k) = g * t0 / gamma * (1 - (mg%d(i,j,k) / pref)**(gamma * Rd / g))
        end do
      end do
    end do
    call fill_halo(gz)

    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gz_lev%d(i,j,k) = g * t0 / gamma * (1 - (mg_lev%d(i,j,k) / pref)**(gamma * Rd / g))
        end do
      end do
    end do
    call fill_halo(gz_lev)
    end associate

  end subroutine rossby_haurwitz_wave_3d_test_set_ic

end module rossby_haurwitz_wave_3d_test_mod
