module albedo_mod

  use const_mod

  implicit none

  private

  public albedo_ocnice

  real(r8), parameter :: ice_alb_sw = 0.7_r8 ! Sea ice albedo for shortwave 0.2-0.7 micro-meters
  real(r8), parameter :: ice_alb_lw = 0.5_r8 ! Sea ice albedo for longwave 0.7-5.0 micro-meters

contains

  subroutine albedo_ocnice(ncol, lndfrac, ocnfrac, icefrac, coszrs, asdir, asdif, aldir, aldif)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: lndfrac(ncol)
    real(r8), intent(in) :: ocnfrac(ncol)
    real(r8), intent(in) :: icefrac(ncol)
    real(r8), intent(in) :: coszrs(ncol)
    real(r8), intent(out) :: asdir(ncol)
    real(r8), intent(out) :: asdif(ncol)
    real(r8), intent(out) :: aldir(ncol)
    real(r8), intent(out) :: aldif(ncol)

    integer i

    do i = 1, ncol
      ! Initialize all ocean/sea ice albedos to 1.
      if (ocnfrac(i) > 0 .or. icefrac(i) > 0) then
        asdir(i) = 1
        asdif(i) = 1
        aldir(i) = 1
        aldif(i) = 1
      end if
      ! Ice-free ocean albedo is a function of solar zenith angle only.
      if (icefrac(i) == 0 .and. coszrs(i) > 0) then
        aldir(i) = (0.026_r8 / (coszrs(i)**1.7_r8 + 0.065_r8)) &
                 + (0.15_r8 * (coszrs(i) - 0.1_r8) &
                            * (coszrs(i) - 0.5_r8) &
                            * (coszrs(i) - 1.0_r8))
        asdir(i) = aldir(i)
        aldif(i) = 0.06_r8
        asdif(i) = 0.06_r8
      end if
      if (icefrac(i) == 1 .and. coszrs(i) > 0) then
        asdir(i) = ice_alb_sw
        aldir(i) = ice_alb_lw
        asdif(i) = ice_alb_sw
        aldif(i) = ice_alb_lw
      end if
    end do

  end subroutine albedo_ocnice

end module albedo_mod