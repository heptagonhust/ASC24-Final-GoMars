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

  real(r8), parameter :: z0_pbl = 0.01_r8 ! FIXME: Put it in other place.

contains

  subroutine mars_sfc_nasa_run(ncol, nlev, wsb, t_sfc, z0, pt, z, ustar, ptstar)

    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    real(r8), intent(in ), dimension(ncol     ) :: wsb
    real(r8), intent(in ), dimension(ncol     ) :: t_sfc
    real(r8), intent(in ), dimension(ncol     ) :: z0
    real(r8), intent(in ), dimension(ncol,nlev) :: pt
    real(r8), intent(in ), dimension(ncol,nlev) :: z
    real(r8), intent(out), dimension(ncol     ) :: ustar
    real(r8), intent(out), dimension(ncol     ) :: ptstar

    real(r8) rib, fh, fm, cdh, cdm, rlnzz
    integer icol

    do icol = 1, ncol
      rib = (g * z(icol,nlev) / (pt(icol,nlev)*wsb(icol)**2 + 1.0e-09_r8)) * (pt(icol,nlev) - t_sfc(icol))
      if (rib >= 0) then
        fh = 1.0_r8 / (1 + (15 * rib / sqrt(1 + 5 * rib)))
        fm = 1.0_r8 / (1 + (10 * rib / sqrt(1 + 5 * rib)))
      else
        fh = sqrt(1 - 64 * rib)
        fm = sqrt(1 - 16 * rib)
      end if
      rlnzz = log(z(icol,nlev) / z0(icol))
      cdh = sqrt(fh) * ka / rlnzz
      cdm = fm * (ka / rlnzz)**2
      ustar(icol) = sqrt(cdm) * wsb(icol)
      ptstar(icol) = cdh * (pt(icol,nlev) - t_sfc(icol))
    end do

  end subroutine mars_sfc_nasa_run

end module mars_sfc_nasa_mod
