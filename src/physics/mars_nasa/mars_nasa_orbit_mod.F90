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

module mars_nasa_orbit_mod

  use const_mod

  implicit none

  private

  public mars_nasa_orbit_init
  public solar_dist
  public solar_decl

  ! Convert Mkm to AU
  real(r8), parameter :: au       = 1.0_r8 / 149.597927_r8
  ! Aphelion Sun-Mars distance (Mkm)
  real(r8), parameter :: aphe     = 249.22_r8
  ! Perihelion Sun-Mars distance (Mkm)
  real(r8), parameter :: peri     = 206.66_r8
  ! Obliquity (deg)
  real(r8), parameter :: obliq    = 25.19_r8
  ! Eccentricity (rad)
  real(r8), parameter :: eccen    = (aphe - peri) / (aphe + peri)
  ! Semimajor axis (AU)
  real(r8), parameter :: semia    = 0.5_r8 * (aphe + peri) * (1 - eccen**2) * au
  ! Sol day per Martian year
  real(r8), parameter :: year_sol = 669
  ! Sol day at perihelion
  real(r8), parameter :: peri_sol = 485
  ! Difference of solar longiutde between Ls~0 and perihelion (rad)
  real(r8) peri_dls

contains

  subroutine mars_nasa_orbit_init()

    real(r8) A, E, M, dE

    ! Calculate eccentric anomaly of Ls~0 using Kepler's equation.
    A = (year_sol - peri_sol) / year_sol
    A = pi2 * (A - int(A))
    M = abs(A)
    E = M + eccen * sin(M)
    dE = 1
    do while (abs(dE) > 1.0e-14)
      dE = - (E - eccen * sin(E) - M) / (1 - eccen * cos(E))
      E = E + dE
    end do
    if (A < 0) E = -E

    peri_dls = 2 * atan(sqrt((1 + eccen) / (1 - eccen)) * tan(E / 2.0_r8)) ! Gauss' equation

  end subroutine mars_nasa_orbit_init

  ! Caluclate the distance between Sun and Mars in AU units.
  pure real(r8) function solar_dist(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    real(r8) alpha, beta

    alpha = ls + peri_dls
    beta = 2 * atan(sqrt((1 - eccen) / (1 + eccen)) * tan(alpha / 2.0_r8))

    res = semia * (1 - eccen * cos(beta))

  end function solar_dist

  pure real(r8) function solar_decl(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = asin(sin(ls) * sin(obliq * rad))

  end function solar_decl

end module mars_nasa_orbit_mod
