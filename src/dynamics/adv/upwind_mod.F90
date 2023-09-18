module upwind_mod

  use const_mod

  implicit none

  private

  public upwind1
  public upwind3
  public upwind5

contains

  pure real(r8) function upwind1(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(0:1)
    real(r8), intent(in) :: wgt

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8

    res = c11 * (f(1) + f(0)) + c12 * (f(1) - f(0)) * wgt * dir

  end function upwind1

  pure real(r8) function upwind3(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: wgt
    real(r8), intent(in) :: f(-1:2)

    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8

    res = c31 * (f(1) + f( 0)) &
        + c32 * (f(2) + f(-1)) &
        + c33 * (f(2) - f(-1) - 3 * (f(1) - f(0))) * wgt * dir

  end function upwind3

  pure real(r8) function upwind5(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: wgt
    real(r8), intent(in) :: f(-2:3)

    real(r8), parameter :: c41 = 37.0_r8 / 60.0_r8
    real(r8), parameter :: c42 = -2.0_r8 / 15.0_r8
    real(r8), parameter :: c43 =  1.0_r8 / 60.0_r8
    real(r8), parameter :: c44 = -1.0_r8 / 60.0_r8

    res = c41 * (f(1) + f( 0)) &
        + c42 * (f(2) + f(-1)) &
        + c43 * (f(3) + f(-2)) &
        + c44 * (f(3) - f(-2) - 5 * (f(2) - f(-1)) + 10 * (f(1) - f(0))) * wgt * dir

  end function upwind5

end module upwind_mod
