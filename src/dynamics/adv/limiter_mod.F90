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
    pure real(r8) function slope_interface(fm1, f, fp1)
      import r8
      real(r8), intent(in) :: fm1, f, fp1
    end function slope_interface
  end interface

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

  pure real(r8) function slope_simple(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    res = (fp1 - fm1) * 0.5_r8

  end function slope_simple

  pure real(r8) function slope_mono(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df, df_min, df_max

    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    df_min = 2 * (f - min(fm1, f, fp1))
    df_max = 2 * (max(fm1, f, fp1) - f)
    res = sign(min(abs(df), df_min, df_max), df)

  end function slope_mono

  pure real(r8) function slope_pd(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df

    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    res = sign(min(abs(df), 2 * f), df)

  end function slope_pd

end module limiter_mod
