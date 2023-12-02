module mountain_wave_test_mod

  use flogger
  use namelist_mod
  use const_mod
  use latlon_parallel_mod
  use block_mod
  use formula_mod
  use operators_mod

  implicit none

  private

  public mountain_wave_test_set_ic

  real(r8), parameter :: T0   = 288.d0      ! K
  real(r8), parameter :: h0   = 2000.d0     ! m
  real(r8), parameter :: d    = 1.5e6
  real(r8), parameter :: u0   = 20.d0       ! m s-1
  real(r8), parameter :: lonc = pi05
  real(r8), parameter :: latc = pi / 6.0
  real(r8), parameter :: kap  = 2.d0 / 7.d0
  real(r8), parameter :: psp  = 93000.d0    ! Pa
  real(r8), parameter :: N    = 0.0182      ! s-1

contains

  subroutine mountain_wave_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               mg     => block%dstate(1)%mg    , &
               pt     => block%dstate(1)%pt    , &
               gzs    => block%static%gzs)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          u%d(i,j,k) = u0 * cos_lat
        end do
      end do
    end do
    call fill_halo(u)

    v%d = 0.0

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        full_lon = mesh%full_lon(i)
        r = radius * acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        gzs%d(i,j) = g * h0 * exp(-(r / d)**2)
      end do
    end do
    call fill_halo(gzs)

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        mgs%d(i,j) = psp * exp(-0.5_r8 * radius * N**2 * u0 / g**2 / kap * (u0 / radius + 2.0_r8 * omega) * &
                     (sin_lat**2 - 1.0_r8) - N**2 / g**2 / kap * gzs%d(i,j))
      end do
    end do
    call fill_halo(mgs)

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = modified_potential_temperature(288.0_r8, mg%d(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(pt)

    if (nonhydrostatic) then
      ! FIXME: Calculate tv.
      call calc_gz_lev(block, block%dstate(1))
    end if
    end associate

  end subroutine mountain_wave_test_set_ic

end module mountain_wave_test_mod
