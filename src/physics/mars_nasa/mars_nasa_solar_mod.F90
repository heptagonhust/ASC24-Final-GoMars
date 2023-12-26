! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_solar_mod

  use const_mod
  use mars_nasa_spectra_mod
  use mars_nasa_objects_mod
  use mars_orbit_mod

  implicit none

  private

  public mars_nasa_solar_init
  public mars_nasa_solar_final
  public calc_solar_flux_toa

  real(r8), allocatable, dimension(:) :: sol_flx_1au_spec ! Solar flux within each visible spectral interval at 1AU (W m-2)
  real(r8), allocatable, dimension(:) :: sol_flx_mars_spec
  real(r8) sol_flx_1au
  real(r8) sol_flx_mars

contains

  subroutine mars_nasa_solar_init()

    if (spec_vis%n /= 7) then
      stop 'mars_nasa_solar_mod only matches with visible spectra with 7 bands!'
    end if

    allocate(sol_flx_1au_spec (spec_vis%n))
    allocate(sol_flx_mars_spec(spec_vis%n))

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    sol_flx_1au_spec = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]
    sol_flx_1au = sum(sol_flx_1au_spec)

  end subroutine mars_nasa_solar_init

  subroutine mars_nasa_solar_final()

    if (allocated(sol_flx_1au_spec )) deallocate(sol_flx_1au_spec )
    if (allocated(sol_flx_mars_spec)) deallocate(sol_flx_mars_spec)

  end subroutine mars_nasa_solar_final

  subroutine calc_solar_flux_toa(ls)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    real(r8) d2

    d2 = solar_dist(ls)**2
    sol_flx_mars_spec = sol_flx_1au_spec / d2
    sol_flx_mars = sol_flx_1au / d2

  end subroutine calc_solar_flux_toa

end module mars_nasa_solar_mod
