module aquaplanet_test_mod

  use const_mod

  implicit none

  private

  public aquaplanet_test_set_bc

contains

  elemental subroutine aquaplanet_test_set_bc(lon, lat, sst, lndfrac, ocnfrac, icefrac)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(out) :: sst
    real(r8), intent(out) :: lndfrac
    real(r8), intent(out) :: ocnfrac
    real(r8), intent(out) :: icefrac

    real(r8), parameter :: max_lat = pi / 3.0_r8
    real(r8), parameter :: max_t   = 27.0_r8
    real(r8), parameter :: min_t   = 0
    real(r8) sin_lat

    sin_lat = sin(pi05 * lat / max_lat)

    ! QOBS
    if (abs(lat) < pi / 3.0_r8) then
      sst = 0.5_r8 * (2 - sin_lat**4 - sin_lat**2) * (max_t - min_t) + min_t
    else
      sst = 0
    end if
    sst = sst + 273.15_r8
    lndfrac = 0
    ocnfrac = 1
    icefrac = 0

  end subroutine aquaplanet_test_set_bc

end module aquaplanet_test_mod