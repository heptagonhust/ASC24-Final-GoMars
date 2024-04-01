! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module limiter_mod

  use flogger
  use const_mod
  use namelist_mod
  use process_mod

  implicit none

  private

  public limiter_init
  public slope

  interface
    real(r8) function slope_interface(fm1, f, fp1)
      import r8
      real(r8), intent(in) :: fm1, f, fp1
      !$omp declare target
    end function slope_interface
  end interface
  !$omp declare target(slope_interface)
  procedure(slope_interface), pointer :: slope => null()

contains

  subroutine limiter_init()

    select case (limiter_type)
    case ('none')
      slope => slope_simple
    case ('mono')
      slope => slope_mono
    case ('pd')
      slope => slope_pd
    case default
      call log_error('Invalid limiter_type ' // trim(limiter_type) // '!', pid=proc%id)
    end select

  end subroutine limiter_init

  real(r8) function slope_simple(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1
    !$omp declare target
    res = (fp1 - fm1) * 0.5_r8

  end function slope_simple

  real(r8) function slope_mono(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df, df_min, df_max
    !$omp declare target
    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    df_min = 2 * (f - min(fm1, f, fp1))
    df_max = 2 * (max(fm1, f, fp1) - f)
    res = sign(min(abs(df), df_min, df_max), df)

  end function slope_mono

  real(r8) function slope_pd(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df
    !$omp declare target
    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    res = sign(min(abs(df), 2 * f), df)

  end function slope_pd

end module limiter_mod
