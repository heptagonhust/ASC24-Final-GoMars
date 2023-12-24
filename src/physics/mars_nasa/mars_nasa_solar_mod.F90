module mars_nasa_solar_mod

  use const_mod
  use mars_nasa_spectra_mod
  use mars_nasa_objects_mod

  implicit none

  private

  public mars_nasa_solar_init
  public mars_nasa_solar_final
  public direct_solar_flux

  real(r8), allocatable, dimension(:) :: sol_flx ! Solar flux within each visible spectral interval at 1AU (W m-2)

contains

  subroutine mars_nasa_solar_init()

    if (spec_vis%n /= 7) then
      stop 'mars_nasa_solar_mod only matches with visible spectra with 7 bands!'
    end if

    allocate(sol_flx(spec_vis%n))

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    sol_flx = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]

  end subroutine mars_nasa_solar_init

  subroutine mars_nasa_solar_final()

    if (allocated(sol_flx)) deallocate(sol_flx)

  end subroutine mars_nasa_solar_final

  pure real(r8) function direct_solar_flux() result(res)

  end function direct_solar_flux

end module mars_nasa_solar_mod
