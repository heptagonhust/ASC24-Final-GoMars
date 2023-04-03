module steady_state_test_mod

  use const_mod
  use formula_mod
  use parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public steady_state_test_set_ic

  real(r8), parameter :: alpha = 0.0_r8
  real(r8), parameter :: u0    = 35        ! m s-1
  real(r8), parameter :: t0    = 288       ! K
  real(r8), parameter :: gamma = 0.005     ! K m-1
  real(r8), parameter :: dt    = 4.8e5     ! K
  real(r8), parameter :: eta0  = 0.252
  real(r8), parameter :: etat  = 0.2       ! Tropopause level

contains

  subroutine steady_state_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) etav, eta, tbar, gzbar, sin_lat, cos_lat
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               mg     => block%dstate(1)%mg    , &
               gz_lev => block%dstate(1)%gz_lev, &
               gz     => block%dstate(1)%gz    , &
               t      => block%dstate(1)%t     , &
               pt     => block%dstate(1)%pt    , &
               gzs    => block%static%gzs)
    mgs = 1.0e5_r8
    v   = 0

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          u(i,j,k) = u0 * cos(etav)**(1.5d0) * sin(2 * mesh%full_lat(j))**2
        end do
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        tbar = t0 * eta**(Rd * gamma / g)
      else
        tbar = t0 * eta**(Rd * gamma / g) + dt * (etat - eta)**5
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          t(i,j,k) = tbar + 3d0 / 4d0 * eta * pi * u0 / Rd * sin(etav) * sqrt(cos(etav)) * (         &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * 2 * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega          &
          )
          pt(i,j,k) = potential_temperature(t(i,j,k), mg(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(block%halo, t, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%half_kds, mesh%half_kde
      eta = mesh%half_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (   &
            (log(eta / etat) + 137d0 / 60d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10d0 / 3d0 * etat**2 * eta**3 +           &
            5d0 / 4d0 * etat * eta**4 - 1d0 / 5d0 * eta**5                   &
          )
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          gz_lev(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                      &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
          )
        end do
      end do
    end do
    call fill_halo(block%halo, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    do k = mesh%full_kds, mesh%full_kde
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (   &
            (log(eta / etat) + 137d0 / 60d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10d0 / 3d0 * etat**2 * eta**3 +           &
            5d0 / 4d0 * etat * eta**4 - 1d0 / 5d0 * eta**5                   &
          )
      end if
      do j = mesh%full_jds, mesh%full_jde
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          gz(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                          &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
          )
        end do
      end do
    end do
    call fill_halo(block%halo, gz, full_lon=.true., full_lat=.true., full_lev=.true.)

    etav = (1 - eta0) * pi / 2
    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        gzs(i,j) = u0 * cos(etav)**1.5d0 * (                                                   &
          (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
          (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
        )
      end do
    end do
    call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine steady_state_test_set_ic

end module steady_state_test_mod
