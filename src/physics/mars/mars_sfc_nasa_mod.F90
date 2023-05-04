! Reference:
!
! - Savijarvi, H. (1995): Mars Boundary Layer Modeling: Diurnal Moisture Cycle
!   and Soil Properties at the Viking Lander 1 Site.
! - Hourdin et al., (1995): The sensitivity of the Martian surface pressure and
!   atmospheric mass budget to various parameters: A comparison between numerical
!   simulations and Viking observations.

module mars_sfc_nasa_mod

  use const_mod

  implicit none

  private

  public mars_sfc_nasa_run

  real(r8), parameter :: z0_pbl = 0.01_r8 ! FIXME: Put it in other place.

contains

  subroutine mars_sfc_nasa_run(ncol, nlev, wsb, t_sfc, z0, pt, z, ustar, ptstar)

    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    real(r8), intent(in ), dimension(ncol     ) :: wsb    ! Wind speed of lowest full level (m/s)
    real(r8), intent(in ), dimension(ncol     ) :: t_sfc  ! Surface temperature (K)
    real(r8), intent(in ), dimension(ncol     ) :: z0     ! Roughness length (m)
    real(r8), intent(in ), dimension(ncol,nlev) :: pt     ! Potential temperature (K)
    real(r8), intent(in ), dimension(ncol,nlev) :: z      ! Geopotential height (m)
    real(r8), intent(out), dimension(ncol     ) :: ustar  ! Friction velocity (m/s)
    real(r8), intent(out), dimension(ncol     ) :: ptstar ! Friction temperature (K)

    real(r8) ri, fh, fm, cdh, cdm, lnz
    integer icol

    do icol = 1, ncol
      ri = (g * z(icol,nlev) * (pt(icol,nlev) - t_sfc(icol)) / (pt(icol,nlev) * wsb(icol)**2 + 1.0e-09_r8))
      if (ri >= 0) then
        fh = 1.0_r8 / (1 + (15 * ri / sqrt(1 + 5 * ri)))
        fm = 1.0_r8 / (1 + (10 * ri / sqrt(1 + 5 * ri)))
      else
        fh = sqrt(1 - 64 * ri)
        fm = sqrt(1 - 16 * ri)
      end if
      lnz = log(z(icol,nlev) / z0(icol))
      cdh = sqrt(fh) * ka / lnz
      cdm = fm * (ka / lnz)**2
      ustar(icol) = sqrt(cdm) * wsb(icol)
      ptstar(icol) = cdh * (pt(icol,nlev) - t_sfc(icol))
    end do

  end subroutine mars_sfc_nasa_run

end module mars_sfc_nasa_mod
