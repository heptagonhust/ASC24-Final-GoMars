module mars_nasa_orbit_mod

  use const_mod

  implicit none

  private

  public mars_nasa_orbit_init

  real(r8), parameter :: au_mkm   = 149.597927_r8 ! AU distance in Mkm
  real(r8), parameter :: aphe     = 249.22_r8     ! Aphelion Sun-Mars distance (Mkm)
  real(r8), parameter :: peri     = 206.66_r8     ! Perihelion Sun-Mars distance (Mkm)
  real(r8), parameter :: obliq    = 25.19_r8      ! Obliquity (deg) 
  real(r8), parameter :: eccen    = (aphe - peri) / (aphe + peri)
  real(r8), parameter :: semia    = 0.5_r8 * (aphe + peri) * (1 - eccen**2) / au_mkm
  real(r8), parameter :: year_sol = 669           ! Sol day per Martian year
  real(r8), parameter :: peri_sol = 485           ! Sol day at perihelion
  real(r8) tperi ! One planet moving period between two perihelion.

contains

  subroutine mars_nasa_orbit_init()

    real(r8) A, E, M, dE

    ! Calculate eccentric anomaly.
    A = (year_sol - peri_sol) / year_sol
    A = pi2 * (A - nint(A))
    M = abs(A)
    E = M + eccen * sin(M)
    dE = 1
    do while (abs(dE) > eps)
      dE = - (E - eccen * sin(E) - M) / (1 - eccen * cos(E))
      E = E + dE
    end do
    if (A < 0) E = -E

    tperi = 2 * atan(sqrt((1 + eccen) / (1 - eccen)) * tan(E / 2.0_r8))

  end subroutine mars_nasa_orbit_init

  pure real(r8) function solar_dist(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = semia / (1 + eccen * cos(ls))

  end function solar_dist

end module mars_nasa_orbit_mod
