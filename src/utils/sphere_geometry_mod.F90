module sphere_geometry_mod

  use flogger
  use const_mod
  use math_mod

  implicit none

  private

  public cartesian_transform
  public inverse_cartesian_transform
  public rotation_transform
  public inverse_rotation_transform
  public calc_distance
  public calc_area
  public calc_area_with_last_small_arc
  public norm_vector
  public calc_sphere_angle
  public calc_arc_length

  integer, parameter :: ORIENT_LEFT  = 1
  integer, parameter :: ORIENT_RIGHT = 2
  integer, parameter :: ORIENT_ON    = 3

  interface cartesian_transform
    module procedure cartesian_transform_1_r4
    module procedure cartesian_transform_2_r4
    module procedure cartesian_transform_1_r8
    module procedure cartesian_transform_2_r8
  end interface cartesian_transform

  interface inverse_cartesian_transform
    module procedure inverse_cartesian_transform_1_r8
    module procedure inverse_cartesian_transform_1_r16
  end interface inverse_cartesian_transform

  interface rotation_transform
    module procedure rotation_transform_r4
    module procedure rotation_transform_r8
    module procedure rotation_transform_r16
  end interface rotation_transform

  interface inverse_rotation_transform
    module procedure inverse_rotation_transform_r4
    module procedure inverse_rotation_transform_r8
    module procedure inverse_rotation_transform_r16
  end interface inverse_rotation_transform

  interface calc_distance
    module procedure calc_distance_r4
    module procedure calc_distance_r8
    module procedure calc_distance_r16
  end interface

  interface calc_sphere_angle
    module procedure calc_sphere_angle_1
  end interface calc_sphere_angle

  interface calc_arc_length
    module procedure calc_arc_length_1
  end interface calc_arc_length

contains

  subroutine cartesian_transform_1_r4(lon, lat, x, y, z)

    real(4), intent(in) :: lon, lat
    real(4), intent(out) :: x, y, z

    real(4) cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_1_r4

  subroutine cartesian_transform_2_r4(lon, lat, x, y, z)

    real(4), intent(in) :: lon, lat
    real(16), intent(out) :: x, y, z

    real(16) cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_2_r4

  subroutine cartesian_transform_1_r8(lon, lat, x, y, z)

    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: x, y, z

    real(8) cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_1_r8

  subroutine cartesian_transform_2_r8(lon, lat, x, y, z)

    real(8), intent(in) :: lon, lat
    real(16), intent(out) :: x, y, z

    real(16) cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_2_r8

  subroutine inverse_cartesian_transform_1_r8(lon, lat, x, y, z)

    real(8), intent(out) :: lon, lat
    real(16), intent(in)  :: x, y, z

    lon = atan2(y, x)
    lat = asin(z / radius)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine inverse_cartesian_transform_1_r8

  subroutine inverse_cartesian_transform_1_r16(lon, lat, x, y, z)

    real(16), intent(out) :: lon, lat
    real(16), intent(in)  :: x, y, z

    lon = atan2(y, x)
    lat = asin(z / radius)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine inverse_cartesian_transform_1_r16

  ! ************************************************************************** !
  ! Rotation transform                                                         !
  ! Purpose:                                                                   !
  !   Calculate the rotating transformation and its inverse of the original    !
  !   coordinate system (lon_o,lat_o) to the rotated one (lon_r, lat_r) with   !
  !   the north pole (lon_p,lat_p) defined at the original coordinate system.  !
  ! ************************************************************************** !

  subroutine rotation_transform_r4(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(4), intent(in) :: lon_p, lat_p ! Rotated pole coordinate
    real(4), intent(in) :: lon_o, lat_o ! Original coordinate
    real(4), intent(out), optional :: lon_r, lat_r ! Rotated coordinate

    real(4) tmp1, tmp2, tmp3, dlon

    dlon = lon_o - lon_p
    if (present(lon_r)) then
        tmp1 = cos(lat_o) * sin(dlon)
        tmp2 = cos(lat_o) * sin(lat_p) * cos(dlon) - cos(lat_p) * sin(lat_o)
        lon_r = atan2(tmp1, tmp2)
        if (lon_r < 0.0d0) lon_r = pi2 + lon_r
    end if
    if (present(lat_r)) then
        tmp1 = sin(lat_o) * sin(lat_p)
        tmp2 = cos(lat_o) * cos(lat_p) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        lat_r = asin(tmp3)
    end if

  end subroutine rotation_transform_r4

  subroutine rotation_transform_r8(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(8), intent(in) :: lon_p, lat_p ! Rotated pole coordinate
    real(8), intent(in) :: lon_o, lat_o ! Original coordinate
    real(8), intent(out), optional :: lon_r, lat_r ! Rotated coordinate

    real(8) tmp1, tmp2, tmp3, dlon

    dlon = lon_o - lon_p
    if (present(lon_r)) then
        tmp1 = cos(lat_o) * sin(dlon)
        tmp2 = cos(lat_o) * sin(lat_p) * cos(dlon) - cos(lat_p) * sin(lat_o)
        lon_r = atan2(tmp1, tmp2)
        if (lon_r < 0.0d0) lon_r = pi2 + lon_r
    end if
    if (present(lat_r)) then
        tmp1 = sin(lat_o) * sin(lat_p)
        tmp2 = cos(lat_o) * cos(lat_p) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        lat_r = asin(tmp3)
    end if

  end subroutine rotation_transform_r8

  subroutine rotation_transform_r16(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(16), intent(in) :: lon_p, lat_p ! Rotated pole coordinate
    real(16), intent(in) :: lon_o, lat_o ! Original coordinate
    real(16), intent(out), optional :: lon_r, lat_r ! Rotated coordinate

    real(16) tmp1, tmp2, tmp3, dlon

    dlon = lon_o - lon_p
    if (present(lon_r)) then
        tmp1 = cos(lat_o) * sin(dlon)
        tmp2 = cos(lat_o) * sin(lat_p) * cos(dlon) - cos(lat_p) * sin(lat_o)
        lon_r = atan2(tmp1, tmp2)
        if (lon_r < 0.0d0) lon_r = pi2 + lon_r
    end if
    if (present(lat_r)) then
        tmp1 = sin(lat_o) * sin(lat_p)
        tmp2 = cos(lat_o) * cos(lat_p) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        lat_r = asin(tmp3)
    end if

  end subroutine rotation_transform_r16

  subroutine inverse_rotation_transform_r4(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(4), intent(in)  :: lon_p, lat_p ! Rotated pole coordinate
    real(4), intent(out) :: lon_o, lat_o ! Original coordinate
    real(4), intent(in)  :: lon_r, lat_r ! Rotated coordinate

    real(4) sin_lon_r, cos_lon_r, sin_lat_r, cos_lat_r, sin_lat_p, cos_lat_p
    real(4) tmp1, tmp2, tmp3

    sin_lon_r = sin(lon_r)
    cos_lon_r = cos(lon_r)
    sin_lat_r = sin(lat_r)
    cos_lat_r = cos(lat_r)
    sin_lat_p = sin(lat_p)
    cos_lat_p = cos(lat_p)

    tmp1 = cos_lat_r * sin_lon_r
    tmp2 = sin_lat_r * cos_lat_p + cos_lat_r * cos_lon_r * sin_lat_p
    ! This trick is due to the inaccuracy of trigonometry calculation.
    if (abs(tmp2) < eps) tmp2 = 0.0d0
    lon_o = atan2(tmp1, tmp2)
    lon_o = lon_p + lon_o
    if (lon_o > pi2) lon_o = lon_o - pi2
    tmp1 = sin_lat_r * sin_lat_p
    tmp2 = cos_lat_r * cos_lat_p * cos_lon_r
    tmp3 = tmp1 - tmp2
    tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
    lat_o = asin(tmp3)

  end subroutine inverse_rotation_transform_r4

  subroutine inverse_rotation_transform_r8(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real(8), intent(in)  :: lon_p, lat_p ! Rotated pole coordinate
    real(8), intent(out) :: lon_o, lat_o ! Original coordinate
    real(8), intent(in)  :: lon_r, lat_r ! Rotated coordinate

    real(8) sin_lon_r, cos_lon_r, sin_lat_r, cos_lat_r, sin_lat_p, cos_lat_p
    real(8) tmp1, tmp2, tmp3

    sin_lon_r = sin(lon_r)
    cos_lon_r = cos(lon_r)
    sin_lat_r = sin(lat_r)
    cos_lat_r = cos(lat_r)
    sin_lat_p = sin(lat_p)
    cos_lat_p = cos(lat_p)

    tmp1 = cos_lat_r * sin_lon_r
    tmp2 = sin_lat_r * cos_lat_p + cos_lat_r * cos_lon_r * sin_lat_p
    ! This trick is due to the inaccuracy of trigonometry calculation.
    if (abs(tmp2) < eps) tmp2 = 0.0d0
    lon_o = atan2(tmp1, tmp2)
    lon_o = lon_p + lon_o
    if (lon_o > pi2) lon_o = lon_o - pi2
    tmp1 = sin_lat_r * sin_lat_p
    tmp2 = cos_lat_r * cos_lat_p * cos_lon_r
    tmp3 = tmp1 - tmp2
    tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
    lat_o = asin(tmp3)

  end subroutine inverse_rotation_transform_r8

  subroutine inverse_rotation_transform_r16(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

      real(16), intent(in)  :: lon_p, lat_p ! Rotated pole coordinate
      real(16), intent(out) :: lon_o, lat_o ! Original coordinate
      real(16), intent(in)  :: lon_r, lat_r ! Rotated coordinate

      real(16) sin_lon_r, cos_lon_r, sin_lat_r, cos_lat_r, sin_lat_p, cos_lat_p
      real(16) tmp1, tmp2, tmp3

      sin_lon_r = sin(lon_r)
      cos_lon_r = cos(lon_r)
      sin_lat_r = sin(lat_r)
      cos_lat_r = cos(lat_r)
      sin_lat_p = sin(lat_p)
      cos_lat_p = cos(lat_p)

      tmp1 = cos_lat_r * sin_lon_r
      tmp2 = sin_lat_r * cos_lat_p + cos_lat_r * cos_lon_r * sin_lat_p
      ! This trick is due to the inaccuracy of trigonometry calculation.
      if (abs(tmp2) < eps) tmp2 = 0.0d0
      lon_o = atan2(tmp1, tmp2)
      lon_o = lon_p + lon_o
      if (lon_o > pi2) lon_o = lon_o - pi2
      tmp1 = sin_lat_r * sin_lat_p
      tmp2 = cos_lat_r * cos_lat_p * cos_lon_r
      tmp3 = tmp1 - tmp2
      tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
      lat_o = asin(tmp3)

  end subroutine inverse_rotation_transform_r16

  pure real(4) function calc_distance_r4(lon1, lat1, lon2, lat2) result(res)

    real(4), intent(in) :: lon1
    real(4), intent(in) :: lat1
    real(4), intent(in) :: lon2
    real(4), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance_r4

  pure real(8) function calc_distance_r8(lon1, lat1, lon2, lat2) result(res)

    real(8), intent(in) :: lon1
    real(8), intent(in) :: lat1
    real(8), intent(in) :: lon2
    real(8), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance_r8

  pure real(16) function calc_distance_r16(lon1, lat1, lon2, lat2) result(res)

    real(16), intent(in) :: lon1
    real(16), intent(in) :: lat1
    real(16), intent(in) :: lon2
    real(16), intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance_r16

  real(16) function calc_area(x, y, z) result(res)

    real(16), intent(in) :: x(:)
    real(16), intent(in) :: y(:)
    real(16), intent(in) :: z(:)

    integer n, im1, i, ip1
    real(16) angle

    n = size(x)
#ifndef NDEBUG
    if (n < 3) then
      call log_error('Spherical polygon number is less than 3!', __FILE__, __LINE__)
    end if
#endif
    res = 0.0
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = calc_sphere_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = radius**2 * (res - (n - 2) * pi)

  end function calc_area

  real(16) function calc_area_with_last_small_arc(x, y, z) result(res)

    real(16), intent(in) :: x(:)
    real(16), intent(in) :: y(:)
    real(16), intent(in) :: z(:)

    integer n
    real(16) xv(3), yv(3), zv(3)
    real(16) lon0, lat0, lon1, lat1, lon2, lat2
    real(16) dlon
    real(16) area1, area2, area3

    if (size(x) /= 3) call log_error('Only support triangle with last edge as small arc!', __FILE__, __LINE__)

    res = calc_area(x, y, z)

    call inverse_cartesian_transform(lon0, lat0, x(1), y(1), z(1))
    call inverse_cartesian_transform(lon1, lat1, x(2), y(2), z(2))
    call inverse_cartesian_transform(lon2, lat2, x(3), y(3), z(3))
    if (lat1 /= lat2) call log_error('Small arc is not valid!', __FILE__, __LINE__)

    if (lat1 == 0.0) then
      ! Small arc is actually a great arc.
      return
    else
      dlon = merge(lon2 - lon1, lon1 - lon2, lat0 > lat1)
      if (dlon < 0.0) dlon = dlon + pi2

      xv(1) = 0.0;  yv(1) = 0.0;
      if (lat0 * lat1 >= 0 .and. abs(lat0) > abs(lat1)) then
        ! Point 0 is at the side with the Pole.
        xv(2) = x(2); yv(2) = y(2); zv(2) = z(2)
        xv(3) = x(3); yv(3) = y(3); zv(3) = z(3)
      else
        ! Point 0 is at the opposite hemisphere.
        xv(2) = x(3); yv(2) = y(3); zv(2) = z(3)
        xv(3) = x(2); yv(3) = y(2); zv(3) = z(2)
      end if
      if (lat1 > 0.0) then
        ! Small arc is at the North Sphere.
        zv(1) = radius
        area1 = radius**2 * dlon * (1.0 - sin(lat1))
      else
        ! Small arc is at the South Sphere.
        zv(1) = -radius
        area1 = radius**2 * dlon * (sin(lat1) + 1.0)
      end if
      area2 = calc_area(xv, yv, zv)
      area3 = area1 - area2
      if (area3 < 0.0 .and. abs(area3) > 1.0e-10) then
        area3 = 0
        !call log_warning('Lune area is negative!', __FILE__, __LINE__)
      end if

      if (lat0 * lat1 >= 0 .and. abs(lat0) > abs(lat1)) then
        res = res + area3
      else
        res = res - area3
      end if
      if (res < 0.0) call log_error('Failed to calculate area with small arc!', __FILE__, __LINE__)
    end if

  end function calc_area_with_last_small_arc

  function norm_vector(x) result(res)

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

  ! Calculate the dihedra angle between plane AB and plane BC.

  real(16) function calc_sphere_angle_1(a, b, c) result(res)

    real(16), intent(in) :: a(3)
    real(16), intent(in) :: b(3)
    real(16), intent(in) :: c(3)

    real(16) nab(3) ! Normal vector of plane AB
    real(16) nbc(3) ! Normal vector of plane BC

    nab = norm_vector(cross_product(a, b))
    nbc = norm_vector(cross_product(b, c))
    res = acos(- max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), a) < 0.0) res = pi2 - res

  end function calc_sphere_angle_1

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  real(16) function calc_arc_length_1(a, b) result(res)

    real(16), intent(in) :: a(3)
    real(16), intent(in) :: b(3)

    res = acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function calc_arc_length_1

end module sphere_geometry_mod
