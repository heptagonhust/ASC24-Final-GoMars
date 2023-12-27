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

module mars_orbit_mod

  use const_mod

  implicit none

  private

  public mars_orbit_init
  public solar_dist
  public solar_decl_angle
  public update_solar_decl_angle
  public solar_cos_zenith_angle

  ! Convert Mkm to AU
  real(r8), parameter :: au        = 1.0_r8 / 149.597927_r8
  ! Aphelion Sun-Mars distance (Mkm)
  real(r8), parameter :: raphe     = 249.22_r8
  ! Perihelion Sun-Mars distance (Mkm)
  real(r8), parameter :: rperi     = 206.66_r8
  ! Obliquity (deg)
  real(r8), parameter :: obliq     = 25.2193_r8 * rad
  real(r8), parameter :: sin_obliq = sin(obliq)
  ! Eccentricity (rad)
  real(r8), parameter :: eccen     = (raphe - rperi) / (raphe + rperi)
  ! Semimajor axis (AU)
  real(r8), parameter :: semia     = rperi / (1 - eccen) * au
  ! Semi-laus rectum
  real(r8), parameter :: pelip     = 0.5_r8 * (raphe + rperi) * (1 - eccen**2) * au
  ! Sol day per Martian year
  real(r8), parameter :: year_sol  = 669
  ! Sol day at perihelion
  real(r8), parameter :: peri_sol  = 485
  ! Difference of solar longiutde between Ls~0 and perihelion (rad)
  real(r8) peri_dls

  real(r8) decl_angle
  real(r8) cos_decl_angle
  real(r8) sin_decl_angle

contains

  subroutine mars_orbit_init()

    real(r8) A, E

    ! Calculate eccentric anomaly of Ls~0 using Kepler's equation.
    A = (year_sol - peri_sol) / year_sol
    A = pi2 * (A - int(A))
    E = eccen_anomaly(abs(A))

    peri_dls = 2 * atan(sqrt((1 + eccen) / (1 - eccen)) * tan(E / 2.0_r8)) ! Gauss' equation

  end subroutine mars_orbit_init

  pure real(r8) function eccen_anomaly(ls) result(E)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    real(r8) dE

    E  = ls + eccen * sin(ls)
    dE = 1
    do while (abs(dE) > 1.0e-12)
      dE = (E - eccen * sin(E) - ls) / (1 - eccen * cos(E))
      E  = E - dE
    end do

  end function eccen_anomaly

  ! Caluclate the distance between Sun and Mars in AU units.
  pure real(r8) function solar_dist(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = semia * (1 - eccen * cos(eccen_anomaly(ls + peri_dls)))

  end function solar_dist

  pure real(r8) function solar_decl_angle(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = asin(sin(ls) * sin_obliq)

  end function solar_decl_angle

  subroutine update_solar_decl_angle(ls)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    decl_angle = solar_decl_angle(ls)
    cos_decl_angle = cos(decl_angle)
    sin_decl_angle = sin(decl_angle)

  end subroutine update_solar_decl_angle

  pure real(r8) function solar_cos_zenith_angle(lon, lat, hour_utc) result(res)

    real(r8), intent(in) :: lon      ! Longitude (rad)
    real(r8), intent(in) :: lat      ! Latitude (rad)
    real(r8), intent(in) :: hour_utc ! Hour

    res = sin(lat) * sin_decl_angle - cos(lat) * cos_decl_angle * cos(pi2 * hour_utc + lon)
    if (res < 1.0e-5_r8) res = 0

  end function solar_cos_zenith_angle

end module mars_orbit_mod
