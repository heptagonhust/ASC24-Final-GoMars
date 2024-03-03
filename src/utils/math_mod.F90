module math_mod

  use flogger
  use const_mod

  private

  public norm_vector
  public cross_product
  public det
  public math_inv_matrix
  public tridiag_thomas
  public exp_two_values
  public swap_two_values
  public round_robin

  interface cross_product
    module procedure cross_product_r8
    module procedure cross_product_r16
  end interface cross_product

  interface det
    module procedure det_r8
    module procedure det_r16
  end interface

  interface exp_two_values
    module procedure exp_two_values_r4
    module procedure exp_two_values_r8
  end interface exp_two_values

contains

  pure function norm_vector(x) result(res)

    real(16), intent(in) :: x(:)
    real(16) res(size(x))

    real(16) n

    n = sqrt(sum(x * x))
    if (n /= 0) then
      res = x / n
    else
      res = x
    end if

  end function norm_vector

  pure function cross_product_r8(x, y) result(res)

    real(8), intent(in) :: x(3)
    real(8), intent(in) :: y(3)
    real(8) res(3)

    res(1) = x(2) * y(3) - x(3) * y(2)
    res(2) = x(3) * y(1) - x(1) * y(3)
    res(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product_r8

  pure function cross_product_r16(x, y) result(res)

    real(16), intent(in) :: x(3)
    real(16), intent(in) :: y(3)
    real(16) res(3)

    res(1) = x(2) * y(3) - x(3) * y(2)
    res(2) = x(3) * y(1) - x(1) * y(3)
    res(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product_r16

  pure recursive function det_r8(x) result(res)

    real(8), intent(in) :: x(:,:)

    integer n, s, i
    real(8) y(size(x,1)-1,size(x,2)-1)

    n = size(x, 1)
    if (n == 1) then
      res = x(1,1)
    else
      res = 0; s = 1
      do i = 1, n
        y(:,:i-1) = x(2:,:i-1)
        y(:,i:) = x(2:,i+1:)
        res = res + s * x(1,i) * det(y)
        s = -s
      end do
    end if

  end function det_r8

  pure recursive function det_r16(x) result(res)

    real(16), intent(in) :: x(:,:)

    integer n, s, i
    real(16) y(size(x,1)-1,size(x,2)-1)

    n = size(x, 1)
    if (n == 1) then
      res = x(1,1)
    else
      res = 0; s = 1
      do i = 1, n
        y(:,:i-1) = x(2:,:i-1)
        y(:,i:) = x(2:,i+1:)
        res = res + s * x(1,i) * det(y)
        s = -s
      end do
    end if

  end function det_r16

  subroutine math_inv_matrix(n, A, B)

    integer , intent(in   ) :: n
    real(r8), intent(inout) :: A(n,n)
    real(r8), intent(  out) :: B(n,n)

    real(r8) C(n,n)
    integer ipiv(n)
    integer i, j, k

    C = 0.0d0
    do i = 1, n
      C(i,i) = 1.0d0
    end do

    call gaussian_elimination(n, A, ipiv)

    do i = 1, n - 1
      do j = i + 1, n
        do k = 1, n
          C(ipiv(j),k) = C(ipiv(j),k) - A(ipiv(j),i) * C(ipiv(i),k)
        end do
      end do
    end do

    do i = 1, n
      B(n,i) = C(ipiv(n),i) / A(ipiv(n),n)
      do j = n - 1, 1, -1
        B(j,i) = C(ipiv(j),i)
        do k = j + 1, n
          B(j,i) = B(j,i) - A(ipiv(j),k) * B(k,i)
        end do
        B(j,i) = B(j,i) / A(ipiv(j),j)
      end do
    end do

  end subroutine math_inv_matrix

  subroutine gaussian_elimination(n, A, ipiv)

    integer , intent(in   ) :: n
    real(r8), intent(inout) :: A(n,n)
    integer , intent(  out) :: ipiv(n)

    real(r8) c1, c(n), pi, pi1, pj
    integer i, j, k, itmp

    do i = 1, n
      ipiv(i) = i
    end do

    do i = 1, n
      c1 = 0.0d0
      do j = 1, n
        c1 = max(c1, abs(A(i,j)))
      end do
      c(i) = c1
    end do

    do j = 1, n - 1
      pi1 = 0.0
      do i = j, n
        pi = abs(A(ipiv(i),j)) / c(ipiv(i))
        if (pi > pi1) then
          pi1 = pi
          k   = i
        end if
      end do

      itmp    = ipiv(j)
      ipiv(j) = ipiv(k)
      ipiv(k) = itmp
      do i = j + 1, n
        pj = A(ipiv(i),j) / A(ipiv(j),j)
        A(ipiv(i),j) = pj
        do k = j + 1, n
          A(ipiv(i),k) = A(ipiv(i),k) - pj * A(ipiv(j),k)
        end do
      end do
    end do

  end subroutine gaussian_elimination

  subroutine tridiag_thomas(a, b, c, d, x)

    real(r8), intent(in ) :: a(:)
    real(r8), intent(in ) :: b(:)
    real(r8), intent(in ) :: c(:)
    real(r8), intent(in ) :: d(:)
    real(r8), intent(out) :: x(:)

    real(r8) gam(size(x)-1), rho(size(x))
    integer n, i

    n = size(x)
    !  _                                                _   _      _     _      _
    ! |  b(1)  c(1)                                      | | x(1  ) |   | d(1  ) |
    ! |  a(2)  b(2)  c(2)                                | | x(2  ) |   | d(2  ) |
    ! |        a(3)  b(3)  c(3)                          | | x(3  ) |   | d(3  ) |
    ! |           ...  ...  ...                          | | ...    | = | ...    |
    ! |                ...  ...  ...                     | | ...    |   | ...    |
    ! |                          a(n-1)  b(n-1)  c(n-1)  | | x(n-1) |   | d(n-1) |
    ! |_                                 a(n  )  b(n  ) _| |_x(n  )_|   |_d(n  )_|
    ! Turn matrix into upper diagonal form
    !  _                                                _   _      _     _      _
    ! |  1     𝛾(1)                                      | | x(1  ) |   | 𝜌(1  ) |
    ! |  0     1     𝛾(2)                                | | x(2  ) |   | 𝜌(2  ) |
    ! |        0     1     𝛾(3)                          | | x(3  ) |   | 𝜌(3  ) |
    ! |              0  ...  ...  ...                    | | ...    | = | ...    |
    ! |                     ...  ...                     | | ...    |   | ...    |
    ! |                                     1     𝛾(n-1) | | x(n-1) |   | 𝜌(n-1) |
    ! |_                                    0     1     _| |_x(n  )_|   |_𝜌(n  )_|
    gam(1) = c(1) / b(1)
    do i = 2, n - 1
      gam(i) = c(i) / (b(i) - a(i) * gam(i-1))
    end do
    rho(1) = d(1) / b(1)
    do i = 2, n
      rho(i) = (d(i) - a(i) * rho(i-1)) / (b(i) - a(i) * gam(i-1))
    end do

    ! Solve the final equations.
    x(n) = rho(n)
    do i = n - 1, 1, -1
      x(i) = rho(i) - gam(i) * x(i+1)
    end do

  end subroutine tridiag_thomas

  pure real(4) function exp_two_values_r4(val1, val2, x0, x1, x) result(res)

    real(4), intent(in) :: val1
    real(4), intent(in) :: val2
    real(4), intent(in) :: x0
    real(4), intent(in) :: x1
    real(4), intent(in) :: x

    real(4) w

    w = exp((x - x0)**2 * log(eps) / (x1 - x0)**2)
    res = w * val1 + (1 - w) * val2

  end function exp_two_values_r4

  pure real(8) function exp_two_values_r8(val1, val2, x0, x1, x) result(res)

    real(8), intent(in) :: val1
    real(8), intent(in) :: val2
    real(8), intent(in) :: x0
    real(8), intent(in) :: x1
    real(8), intent(in) :: x

    real(8) w

    w = exp((x - x0)**2 * log(eps) / (x1 - x0)**2)
    res = w * val1 + (1 - w) * val2

  end function exp_two_values_r8

  subroutine swap_two_values(x, y)

    integer, intent(inout) :: x
    integer, intent(inout) :: y

    integer tmp

    tmp = x
    x = y
    y = tmp

  end subroutine swap_two_values

  subroutine round_robin(dim, coord, num, ibeg, iend)

    integer, intent(in) :: dim
    integer, intent(in) :: coord
    integer, intent(inout) :: num
    integer, intent(out) :: ibeg ! Start from 1.
    integer, intent(out) :: iend ! Start from 1.

    integer res_num, tmp_num, i

    res_num = mod(num, dim)
    ibeg = 1
    do i = 0, coord - 1
      if (res_num /= 0 .and. i < res_num) then
        tmp_num = num / dim + 1
      else
        tmp_num = num / dim
      end if
      ibeg = ibeg + tmp_num
    end do
    if (res_num /= 0 .and. coord < res_num) then
      num = num / dim + 1
    else
      num = num / dim
    end if
    iend = ibeg + num - 1

  end subroutine round_robin

end module math_mod
