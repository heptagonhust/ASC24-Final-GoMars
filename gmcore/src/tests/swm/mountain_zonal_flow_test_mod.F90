module mountain_zonal_flow_test_mod

  use flogger
  use const_mod
  use latlon_parallel_mod
  use block_mod

  implicit none

  private

  public mountain_zonal_flow_test_set_ic

  real(r8), parameter :: alpha = 0.0
  real(r8), parameter :: u0 = 20.0
  real(r8)            :: gz0
  real(r8), parameter :: lon0 = pi * 1.5
  real(r8), parameter :: lat0 = pi / 6.0
  real(r8)            :: gzs0
  real(r8), parameter :: R = pi / 9.0

contains

  subroutine mountain_zonal_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    integer i, j, k

    gz0  = 5960.0d0 * g
    gzs0 = 2000.0d0 * g

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    associate (mesh   => block%mesh           , &
               u      => block%dstate(1)%u_lon, &
               v      => block%dstate(1)%v_lat, &
               gz     => block%dstate(1)%gz   , &
               gzs    => block%static%gzs)
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        dlon = abs(mesh%full_lon(i) - lon0)
        dlon = min(dlon, 2 * pi - dlon)
        d = min(R, sqrt(dlon**2 + (mesh%full_lat(j) - lat0)**2))
        gzs%d(i,j) = gzs0 * (1.0 - d / R)
      end do
    end do
    call fill_halo(gzs)

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_ids, mesh%half_ide
        cos_lon = mesh%half_cos_lon(i)
        u%d(i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do
    call fill_halo(u)

    do j = mesh%half_jds, mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        sin_lon = mesh%full_sin_lon(i)
        v%d(i,j,1) = - u0 * sin_lon * sin_alpha
      end do
    end do
    call fill_halo(v)

    do j = mesh%full_jds, mesh%full_jde
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        cos_lon = mesh%full_cos_lon(i)
        gz%d(i,j,1) = gz0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do
    call fill_halo(gz)
    end associate

  end subroutine mountain_zonal_flow_test_set_ic

end module mountain_zonal_flow_test_mod
