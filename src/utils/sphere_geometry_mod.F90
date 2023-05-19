module sphere_geometry_mod

  use const_mod
  use math_mod

  implicit none

  private

  public lonlat2xyz
  public xyz2lonlat
  public rotate
  public rotate_back
  public great_circle
  public spherical_area
  public spherical_area_with_last_small_arc
  public spherical_triangle_area
  public spherical_rectangle_area
  public spherical_angle
  public sphere_geometry_error_message

  integer, parameter :: ORIENT_LEFT  = 1
  integer, parameter :: ORIENT_RIGHT = 2
  integer, parameter :: ORIENT_ON    = 3

  integer, parameter :: err_msg_len = 256

  ! Error codes:
  integer, parameter :: ERR_1 = 1
  character(err_msg_len), parameter :: ERR_1_MSG = 'Spherical polygon number is less than 3!'
  integer, parameter :: ERR_2 = 2
  character(err_msg_len), parameter :: ERR_2_MSG = 'Only support triangle with last edge as small arc!'
  integer, parameter :: ERR_3 = 3
  character(err_msg_len), parameter :: ERR_3_MSG = 'Small arc is not valid!'
  integer, parameter :: ERR_4 = 4
  character(err_msg_len), parameter :: ERR_4_MSG = 'Failed to calculate area with small arc!'

  interface lonlat2xyz
    module procedure lonlat2xyz_1_r4
    module procedure lonlat2xyz_2_r4
    module procedure lonlat2xyz_1_r8
    module procedure lonlat2xyz_2_r8
  end interface lonlat2xyz

  interface xyz2lonlat
    module procedure xyz2lonlat_1_r8
    module procedure xyz2lonlat_2_r8
    module procedure xyz2lonlat_1_r16
  end interface xyz2lonlat

  interface rotate
    module procedure rotate_r4
    module procedure rotate_r8
    module procedure rotate_r16
  end interface rotate

  interface rotate_back
    module procedure rotate_back_r4
    module procedure rotate_back_r8
    module procedure rotate_back_r16
  end interface rotate_back

  interface spherical_triangle_area
    module procedure spherical_triangle_area_r8
    module procedure spherical_triangle_area_r16
  end interface spherical_triangle_area

  interface spherical_rectangle_area
    module procedure spherical_rectangle_area_r8
    module procedure spherical_rectangle_area_r16
  end interface spherical_rectangle_area

  interface great_circle
    module procedure great_circle_1_r4
    module procedure great_circle_1_r8
    module procedure great_circle_1_r16
    module procedure great_circle_2_r16
  end interface

  interface spherical_angle
    module procedure spherical_angle_1
  end interface spherical_angle

contains

  subroutine lonlat2xyz_1_r4(R, lon, lat, x, y, z)

    real(r8), intent(in) :: R
    real(4), intent(in) :: lon, lat
    real(4), intent(out) :: x, y, z

    real(4) cos_lat

    cos_lat = cos(lat)
    x = R * cos_lat * cos(lon)
    y = R * cos_lat * sin(lon)
    z = R * sin(lat)

  end subroutine lonlat2xyz_1_r4

  subroutine lonlat2xyz_2_r4(R, lon, lat, x, y, z)

    real(r8), intent(in) :: R
    real(4), intent(in) :: lon, lat
    real(16), intent(out) :: x, y, z

    real(16) cos_lat

    cos_lat = cos(lat)
    x = R * cos_lat * cos(lon)
    y = R * cos_lat * sin(lon)
    z = R * sin(lat)

  end subroutine lonlat2xyz_2_r4

  subroutine lonlat2xyz_1_r8(R, lon, lat, x, y, z)

    real(r8), intent(in) :: R
    real(8), intent(in) :: lon, lat
    real(8), intent(out) :: x, y, z

    real(8) cos_lat

    cos_lat = cos(lat)
    x = R * cos_lat * cos(lon)
    y = R * cos_lat * sin(lon)
    z = R * sin(lat)

  end subroutine lonlat2xyz_1_r8

  subroutine lonlat2xyz_2_r8(R, lon, lat, x, y, z)

    real(r8), intent(in) :: R
    real(8), intent(in) :: lon, lat
    real(16), intent(out) :: x, y, z

    real(16) cos_lat

    cos_lat = cos(lat)
    x = R * cos_lat * cos(lon)
    y = R * cos_lat * sin(lon)
    z = R * sin(lat)

  end subroutine lonlat2xyz_2_r8

  subroutine xyz2lonlat_1_r8(R, x, y, z, lon, lat)

    real(r8), intent(in) :: R
    real(8), intent(in)  :: x, y, z
    real(8), intent(out) :: lon, lat

    lon = atan2(y, x)
    lat = asin(z / R)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine xyz2lonlat_1_r8

  subroutine xyz2lonlat_2_r8(R, x, y, z, lon, lat)

    real(r8), intent(in) :: R
    real(16), intent(in)  :: x, y, z
    real(8), intent(out) :: lon, lat

    lon = atan2(y, x)
    lat = asin(z / R)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine xyz2lonlat_2_r8

  subroutine xyz2lonlat_1_r16(R, x, y, z, lon, lat)

    real(r8), intent(in) :: R
    real(16), intent(in)  :: x, y, z
    real(16), intent(out) :: lon, lat

    lon = atan2(y, x)
    lat = asin(z / R)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine xyz2lonlat_1_r16

  ! ************************************************************************** !
  ! Rotation transform                                                         !
  ! Purpose:                                                                   !
  !   Calculate the rotating transformation and its inverse of the original    !
  !   coordinate system (lon_o,lat_o) to the rotated one (lon_r, lat_r) with   !
  !   the north pole (lon_p,lat_p) defined at the original coordinate system.  !
  ! ************************************************************************** !

  subroutine rotate_r4(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_r4

  subroutine rotate_r8(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_r8

  subroutine rotate_r16(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_r16

  subroutine rotate_back_r4(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_back_r4

  subroutine rotate_back_r8(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_back_r8

  subroutine rotate_back_r16(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

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

  end subroutine rotate_back_r16

  pure real(4) function great_circle_1_r4(R, lon1, lat1, lon2, lat2) result(res)

    real(r8), intent(in) :: R
    real(4), intent(in) :: lon1
    real(4), intent(in) :: lat1
    real(4), intent(in) :: lon2
    real(4), intent(in) :: lat2

    res = R * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function great_circle_1_r4

  pure real(8) function great_circle_1_r8(R, lon1, lat1, lon2, lat2) result(res)

    real(r8), intent(in) :: R
    real(8), intent(in) :: lon1
    real(8), intent(in) :: lat1
    real(8), intent(in) :: lon2
    real(8), intent(in) :: lat2

    res = R * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function great_circle_1_r8

  pure real(16) function great_circle_1_r16(R, lon1, lat1, lon2, lat2) result(res)

    real(r8), intent(in) :: R
    real(16), intent(in) :: lon1
    real(16), intent(in) :: lat1
    real(16), intent(in) :: lon2
    real(16), intent(in) :: lat2

    res = R * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function great_circle_1_r16

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  pure real(16) function great_circle_2_r16(R, a, b) result(res)

    real(r8), intent(in) :: R
    real(16), intent(in) :: a(3)
    real(16), intent(in) :: b(3)

    res = R * acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function great_circle_2_r16

  real(16) function spherical_area(R, x, y, z, ierr) result(res)

    real(r8), intent(in) :: R
    real(16), intent(in) :: x(:)
    real(16), intent(in) :: y(:)
    real(16), intent(in) :: z(:)
    integer, intent(out) :: ierr

    integer n, im1, i, ip1
    real(16) angle

    ierr = 0
    n = size(x)
    if (n < 3) then
      ierr = ERR_1
      return
    end if
    res = 0.0
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = spherical_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = R**2 * (res - (n - 2) * pi)

  end function spherical_area

  real(16) function spherical_area_with_last_small_arc(R, x, y, z, ierr) result(res)

    real(r8), intent(in) :: R
    real(16), intent(in) :: x(:)
    real(16), intent(in) :: y(:)
    real(16), intent(in) :: z(:)
    integer, intent(out) :: ierr

    integer n
    real(16) xv(3), yv(3), zv(3)
    real(16) lon0, lat0, lon1, lat1, lon2, lat2
    real(16) dlon
    real(16) area1, area2, area3

    ierr = 0
    if (size(x) /= 3) then
      ierr = ERR_2
      return
    end if

    res = spherical_area(R, x, y, z, ierr)
    if (ierr /= 0) return

    call xyz2lonlat(R, x(1), y(1), z(1), lon0, lat0)
    call xyz2lonlat(R, x(2), y(2), z(2), lon1, lat1)
    call xyz2lonlat(R, x(3), y(3), z(3), lon2, lat2)
    if (lat1 /= lat2) then
      ierr = ERR_3
      return
    end if

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
      area2 = spherical_area(R, xv, yv, zv, ierr)
      if (ierr /= 0) return
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
      if (res < 0.0) then
        ierr = ERR_4
        return
      end if
    end if

  end function spherical_area_with_last_small_arc

  real(8) function spherical_triangle_area_r8(R, p1, p2, p3) result(res)

    real(8), intent(in) :: R
    real(8), intent(in) :: p1(:)
    real(8), intent(in) :: p2(:)
    real(8), intent(in) :: p3(:)

    real(8), dimension(3) :: p31, p12, p23
    real(8) a312, a123, a231

    p31 = cross_product(p3, p1); p31 = p31 / norm2(p31)
    p12 = cross_product(p1, p2); p12 = p12 / norm2(p12)
    p23 = cross_product(p2, p3); p23 = p23 / norm2(p23)
    a312 = acos(-dot_product(p31, p12))
    a123 = acos(-dot_product(p12, p23))
    a231 = acos(-dot_product(p23, p31))
    res = R**2 * (a312 + a123 + a231 - pi)

  end function spherical_triangle_area_r8

  real(16) function spherical_triangle_area_r16(R, p1, p2, p3) result(res)

    real(16), intent(in) :: R
    real(16), intent(in) :: p1(:)
    real(16), intent(in) :: p2(:)
    real(16), intent(in) :: p3(:)

    real(16), dimension(3) :: p31, p12, p23
    real(16) a312, a123, a231

    p31 = cross_product(p3, p1); p31 = p31 / norm2(p31)
    p12 = cross_product(p1, p2); p12 = p12 / norm2(p12)
    p23 = cross_product(p2, p3); p23 = p23 / norm2(p23)
    a312 = acos(-dot_product(p31, p12))
    a123 = acos(-dot_product(p12, p23))
    a231 = acos(-dot_product(p23, p31))
    res = R**2 * (a312 + a123 + a231 - pi)

  end function spherical_triangle_area_r16

  real(8) function spherical_rectangle_area_r8(R, p1, p2, p3, p4) result(res)

    real(8), intent(in) :: R
    real(8), intent(in) :: p1(:)
    real(8), intent(in) :: p2(:)
    real(8), intent(in) :: p3(:)
    real(8), intent(in) :: p4(:)

    real(8), dimension(3) :: p41, p12, p23, p34
    real(8) a412, a123, a234, a341

    p41 = cross_product(p4, p1); p41 = p41 / norm2(p41)
    p12 = cross_product(p1, p2); p12 = p12 / norm2(p12)
    p23 = cross_product(p2, p3); p23 = p23 / norm2(p23)
    p34 = cross_product(p3, p4); p34 = p34 / norm2(p34)
    a412 = acos(-dot_product(p41, p12))
    a123 = acos(-dot_product(p12, p23))
    a234 = acos(-dot_product(p23, p34))
    a341 = acos(-dot_product(p34, p41))
    res = R**2 * (a412 + a123 + a234 + a341 - pi2)

  end function spherical_rectangle_area_r8

  real(16) function spherical_rectangle_area_r16(R, p1, p2, p3, p4) result(res)

    real(16), intent(in) :: R
    real(16), intent(in) :: p1(:)
    real(16), intent(in) :: p2(:)
    real(16), intent(in) :: p3(:)
    real(16), intent(in) :: p4(:)

    real(16), dimension(3) :: p41, p12, p23, p34
    real(16) a412, a123, a234, a341

    p41 = cross_product(p4, p1); p41 = p41 / norm2(p41)
    p12 = cross_product(p1, p2); p12 = p12 / norm2(p12)
    p23 = cross_product(p2, p3); p23 = p23 / norm2(p23)
    p34 = cross_product(p3, p4); p34 = p34 / norm2(p34)
    a412 = acos(-dot_product(p41, p12))
    a123 = acos(-dot_product(p12, p23))
    a234 = acos(-dot_product(p23, p34))
    a341 = acos(-dot_product(p34, p41))
    res = R**2 * (a412 + a123 + a234 + a341 - pi2)

  end function spherical_rectangle_area_r16

  ! Calculate the dihedra angle between plane AB and plane BC.

  pure real(16) function spherical_angle_1(a, b, c) result(res)

    real(16), intent(in) :: a(3)
    real(16), intent(in) :: b(3)
    real(16), intent(in) :: c(3)

    real(16) nab(3) ! Normal vector of plane AB
    real(16) nbc(3) ! Normal vector of plane BC

    nab = norm_vector(cross_product(a, b))
    nbc = norm_vector(cross_product(b, c))
    res = acos(-max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), a) < 0.0) res = pi2 - res

  end function spherical_angle_1

  character(err_msg_len) function sphere_geometry_error_message(ierr) result(res)

    integer, intent(in) :: ierr

    select case (ierr)
    case (0)
      res = 'No error.'
    case (1)
      res = ERR_1_MSG
    case (2)
      res = ERR_2_MSG
    case (3)
      res = ERR_3_MSG
    case (4)
      res = ERR_4_MSG
    case default
      res = 'Unknown error code!'
    end select

  end function sphere_geometry_error_message

end module sphere_geometry_mod
