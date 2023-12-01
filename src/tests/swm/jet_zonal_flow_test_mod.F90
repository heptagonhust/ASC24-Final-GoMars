module jet_zonal_flow_test_mod

  use flogger
  use string
  use const_mod
  use latlon_parallel_mod
  use block_mod

  implicit none

  private

  public jet_zonal_flow_test_set_ic

  real(8), parameter :: u_max = 80.0d0
  real(8), parameter :: lat0 = pi / 7.0d0
  real(8), parameter :: lat1 = pi / 2.0d0 - lat0
  real(8), parameter :: en = exp(-4.0d0 / (lat1 - lat0)**2)
  real(8)            :: gh0
  real(8)            :: ghd
  real(8), parameter :: lat2 = pi / 4.0d0
  real(8), parameter :: alpha = 1.0d0 / 3.0d0
  real(8), parameter :: beta = 1.0d0 / 15.0d0

contains

  subroutine jet_zonal_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, neval, ierr
    real(8) gz_, abserr

    associate (mesh   => block%mesh           , &
               u      => block%dstate(1)%u_lon, &
               v      => block%dstate(1)%v_lat, &
               gz     => block%dstate(1)%gz   , &
               gzs    => block%static%gzs)
    gh0 = g * 1.0d4
    ghd = g * 120.0d0

    gzs%d = 0

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%half_ids, mesh%half_ide
        u%d(i,j,1) = u_function(mesh%full_lat(j))
      end do
    end do
    call fill_halo(u)

    v%d = 0

    do j = mesh%full_jds, mesh%full_jde
      i = mesh%half_ids
      if (j == 1) then
        gz%d(i,j,1) = gh0
      else
        call qags(gh_integrand, -0.5d0*pi, mesh%full_lat(j), 1.0d-12, 1.0d-3, gz_, abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_str(i) // ',' // to_str(j) // ')!', __FILE__, __LINE__)
        end if
        gz%d(i,j,1) = gh0 - gz_
      end if
      do i = mesh%half_ids, mesh%half_ide
        gz%d(i,j,1) = gz%d(mesh%half_ids,j,1)
        ! Add perturbation.
        gz%d(i,j,1) = gz%d(i,j,1) + ghd * &
          cos(mesh%full_lat(j)) * &
          exp(-(merge(mesh%full_lon(i) - 2*pi, mesh%full_lon(i), mesh%full_lon(i) > pi)  / alpha)**2) * &
          exp(-((lat2 - mesh%full_lat(j)) / beta)**2)
      end do
    end do
    call fill_halo(gz)
    end associate

  end subroutine jet_zonal_flow_test_set_ic

  real(8) function gh_integrand(lat) result(res)

    real(8), intent(in) :: lat

    real(8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(8) function u_function(lat) result(res)

    real(8), intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
