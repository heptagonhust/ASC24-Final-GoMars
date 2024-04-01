! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module ppm_mod

  use const_mod
  use namelist_mod
  ! use limiter_mod

  implicit none

  private

  public ppm

contains
  real(r8) function slope(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df, df_min, df_max
    !$omp declare target
    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    df_min = 2 * (f - min(fm1, f, fp1))
    df_max = 2 * (max(fm1, f, fp1) - f)
    res = sign(min(abs(df), df_min, df_max), df)

  end function slope

  subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)

    real(r8), intent(in ) :: fm2
    real(r8), intent(in ) :: fm1
    real(r8), intent(in ) :: f
    real(r8), intent(in ) :: fp1
    real(r8), intent(in ) :: fp2
    real(r8), intent(out) :: fl
    real(r8), intent(out) :: df
    real(r8), intent(out) :: f6
    !$omp declare target

    real(r8) dfl, dfr, fr
    ! interface
    !   real(r8) function slope(fm1, f, fp1)
    !     import r8
    !     real(r8), intent(in) :: fm1, f, fp1
    !     !$omp declare target
    !   end function slope
    ! end interface

    ! Calculate values at left and right cell interfaces.
    dfl = slope(fm2, fm1, f  )
    df  = slope(fm1, f  , fp1)
    dfr = slope(f  , fp1, fp2)
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5_r8 * (fm1 + f) + (dfl - df) / 6.0_r8
    fr = 0.5_r8 * (fp1 + f) + (df - dfr) / 6.0_r8
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f - sign(min(abs(df), abs(fl - f)), df)
    fr = f + sign(min(abs(df), abs(fr - f)), df)
    f6 = 6 * f - 3 * (fl + fr)
    df = fr - fl

  end subroutine ppm

end module ppm_mod
