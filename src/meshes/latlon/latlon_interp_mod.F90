! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_interp_mod

  use flogger
  use const_mod
  use latlon_mesh_mod
  use process_mod
  use latlon_parallel_mod

  implicit none

  private

  public latlon_interp_bilinear_cell
  public latlon_interp_bilinear_lon_edge
  public latlon_interp_bilinear_lat_edge
  public latlon_interp_bilinear_column

contains

  subroutine latlon_interp_bilinear_cell(src_lon, src_lat, src_data, dst_mesh, dst_data, extrap, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(latlon_mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%full_ims:dst_mesh%full_ime, &
                                      dst_mesh%full_jms:dst_mesh%full_jme)
    logical, intent(in), optional :: extrap
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2, pole
    logical is_found, extrap_opt, zero_pole_opt

    extrap_opt = .true.; if (present(extrap)) extrap_opt = extrap
    zero_pole_opt = .false.; if (present(zero_pole)) zero_pole_opt = zero_pole

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%full_ids
    ub_x2 = dst_mesh%full_ide
    lb_y2 = dst_mesh%full_jds
    ub_y2 = dst_mesh%full_jde

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%full_lon_deg, &
               y2 => dst_mesh%full_lat_deg)
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!', __FILE__, __LINE__)
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if (y2(j) < y1(jj)) exit
          if (y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (extrap_opt) then
            if (y2(j) < y1(jj)) then ! south
              j1(j) = jj
              j2(j) = jj + 1
              tmp1 = y1(j1(j)) - y2(j)
              tmp2 = y1(j2(j)) - y2(j)
            else ! north
              j1(j) = jj
              j2(j) = jj - 1
              tmp1 = y2(j) - y1(j1(j))
              tmp2 = y2(j) - y1(j2(j))
            end if
            ywgt1(j) =  tmp2 / (tmp2 - tmp1)
            ywgt2(j) = -tmp1 / (tmp2 - tmp1)
          else if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!', __FILE__, __LINE__)
          end if
        else
          tmp1 = y1(j1(j))
          tmp2 = y1(j2(j))
          ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
          ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
        end if
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%full_jds) = 0.0_r8
      else
        call zonal_sum(proc%zonal_circle, dst_data(dst_mesh%full_ids:dst_mesh%full_ide,dst_mesh%full_jds), pole)
        dst_data(:,dst_mesh%full_jds) = pole / global_mesh%full_nlon
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%full_jde) = 0.0_r8
      else
        call zonal_sum(proc%zonal_circle, dst_data(dst_mesh%full_ids:dst_mesh%full_ide,dst_mesh%full_jde), pole)
        dst_data(:,dst_mesh%full_jde) = pole / global_mesh%full_nlon
      end if
    end if

  end subroutine latlon_interp_bilinear_cell

  subroutine latlon_interp_bilinear_lon_edge(src_lon, src_lat, src_data, dst_mesh, dst_data, extrap, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(latlon_mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%half_ims:dst_mesh%half_ime, &
                                      dst_mesh%full_jms:dst_mesh%full_jme)
    logical, intent(in), optional :: extrap
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2, pole
    logical is_found, extrap_opt, zero_pole_opt

    extrap_opt = .true.; if (present(extrap)) extrap_opt = extrap
    zero_pole_opt = .false.; if (present(zero_pole)) zero_pole_opt = zero_pole

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%half_ids
    ub_x2 = dst_mesh%half_ide
    lb_y2 = dst_mesh%full_jds
    ub_y2 = dst_mesh%full_jde

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%half_lon_deg, &
               y2 => dst_mesh%full_lat_deg)
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!', __FILE__, __LINE__)
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if (y2(j) < y1(jj)) exit
          if (y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (extrap_opt) then
            if (y2(j) < y1(jj)) then ! south
              j1(j) = jj
              j2(j) = jj + 1
              tmp1 = y1(j1(j)) - y2(j)
              tmp2 = y1(j2(j)) - y2(j)
            else ! north
              j1(j) = jj
              j2(j) = jj - 1
              tmp1 = y2(j) - y1(j1(j))
              tmp2 = y2(j) - y1(j2(j))
            end if
            ywgt1(j) =  tmp2 / (tmp2 - tmp1)
            ywgt2(j) = -tmp1 / (tmp2 - tmp1)
          else if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!')
          end if
        else
          tmp1 = y1(j1(j))
          tmp2 = y1(j2(j))
          ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
          ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
        end if
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%full_jds) = 0.0_r8
      else
        call zonal_sum(proc%zonal_circle, dst_data(dst_mesh%half_ids:dst_mesh%half_ide,dst_mesh%full_jds), pole)
        dst_data(:,dst_mesh%full_jds) = pole / global_mesh%half_nlon
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%full_jde) = 0.0_r8
      else
        call zonal_sum(proc%zonal_circle, dst_data(dst_mesh%half_ids:dst_mesh%half_ide,dst_mesh%full_jde), pole)
        dst_data(:,dst_mesh%full_jde) = pole / global_mesh%half_nlon
      end if
    end if

  end subroutine latlon_interp_bilinear_lon_edge

  subroutine latlon_interp_bilinear_lat_edge(src_lon, src_lat, src_data, dst_mesh, dst_data, extrap, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(latlon_mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%full_ims:dst_mesh%full_ime, &
                                      dst_mesh%half_jms:dst_mesh%half_jme)
    logical, intent(in), optional :: extrap
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2
    logical is_found, extrap_opt, zero_pole_opt

    extrap_opt = .true.; if (present(extrap)) extrap_opt = extrap
    zero_pole_opt = .false.; if (present(zero_pole)) zero_pole_opt = zero_pole

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%full_ids
    ub_x2 = dst_mesh%full_ide
    lb_y2 = dst_mesh%half_jds
    ub_y2 = dst_mesh%half_jde

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%full_lon_deg, &
               y2 => dst_mesh%half_lat_deg)
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!')
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if (y2(j) < y1(jj)) exit
          if (y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (extrap_opt) then
            if (y2(j) < y1(jj)) then ! south
              j1(j) = jj
              j2(j) = jj + 1
              tmp1 = y1(j1(j)) - y2(j)
              tmp2 = y1(j2(j)) - y2(j)
            else ! north
              j1(j) = jj
              j2(j) = jj - 1
              tmp1 = y2(j) - y1(j1(j))
              tmp2 = y2(j) - y1(j2(j))
            end if
            ywgt1(j) =  tmp2 / (tmp2 - tmp1)
            ywgt2(j) = -tmp1 / (tmp2 - tmp1)
          else if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!')
          end if
        else
          tmp1 = y1(j1(j))
          tmp2 = y1(j2(j))
          ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
          ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
        end if
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%half_jds) = 0.0_r8
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (zero_pole_opt) then
        dst_data(:,dst_mesh%half_jde) = 0.0_r8
      end if
    end if

  end subroutine latlon_interp_bilinear_lat_edge

  subroutine latlon_interp_bilinear_column(src_lon, src_lat, src_data, dst_lon, dst_lat, dst_data, extrap, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    real(r8), intent(in) :: dst_lon(:)
    real(r8), intent(in) :: dst_lat(:)
    real(r8), intent(out) :: dst_data(:)
    logical, intent(in), optional :: extrap
    integer, intent(out), optional :: ierr

    integer nx1, ny1, nc2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2
    logical is_found, extrap_opt, zero_pole_opt

    extrap_opt = .true.; if (present(extrap)) extrap_opt = extrap

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    nc2 = size(dst_lon)

    allocate(i1(1:nc2), i2(1:nc2))
    allocate(j1(1:nc2), j2(1:nc2))
    allocate(xwgt1(1:nc2), xwgt2(1:nc2))
    allocate(ywgt1(1:nc2), ywgt2(1:nc2))

    associate (x1 => src_lon, &
               y1 => src_lat, &
               x2 => dst_lon, &
               y2 => dst_lat)
    do i = 1, nc2
      is_found = .false.
      do ii = 1, nx1 - 1
        if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
          i1(i) = ii
          i2(i) = ii+1
          tmp1 = x1(i1(i))
          tmp2 = x1(i2(i))
          is_found = .true.
          exit
        else if (x2(i) >= x1(nx1)) then
          i1(i) = nx1
          i2(i) = 1
          tmp1 = x1(i1(i))
          tmp2 = x1(i2(i)) + 360.0d0
          is_found = .true.
          exit
        end if
      end do
      if (.not. is_found) then
        if (present(ierr)) then
          ierr = -1
          deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
          return
        else
          call log_error('latlon_interp_bilinear: Longitude mismatch!', __FILE__, __LINE__)
        end if
      end if
      xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
      xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
    end do

    do j = 1, nc2
      is_found = .false.
      do jj = 1, ny1-1
        if (y2(j) < y1(jj)) exit
        if (y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) then
          j1(j) = jj
          j2(j) = jj + 1
          is_found = .true.
          exit
        end if
      end do
      if (.not. is_found) then
        if (extrap_opt) then
          if (y2(j) < y1(jj)) then ! south
            j1(j) = jj
            j2(j) = jj + 1
            tmp1 = y1(j1(j)) - y2(j)
            tmp2 = y1(j2(j)) - y2(j)
          else ! north
            j1(j) = jj
            j2(j) = jj - 1
            tmp1 = y2(j) - y1(j1(j))
            tmp2 = y2(j) - y1(j2(j))
          end if
          ywgt1(j) =  tmp2 / (tmp2 - tmp1)
          ywgt2(j) = -tmp1 / (tmp2 - tmp1)
        else if (present(ierr)) then
          ierr = -1
          deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
          return
        else
          call log_error('latlon_interp_bilinear: Latitude mismatch!', __FILE__, __LINE__)
        end if
      else
        tmp1 = y1(j1(j))
        tmp2 = y1(j2(j))
        ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
        ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
      end if
    end do

    do i = 1, nc2
      tmp1 = src_data(i1(i),j1(i)) * xwgt1(i) + src_data(i2(i),j1(i)) * xwgt2(i)
      tmp2 = src_data(i1(i),j2(i)) * xwgt1(i) + src_data(i2(i),j2(i)) * xwgt2(i)
      dst_data(i) = tmp1 * ywgt1(i) + tmp2 * ywgt2(i)
    end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

  end subroutine latlon_interp_bilinear_column

end module latlon_interp_mod
