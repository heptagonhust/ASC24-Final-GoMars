! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_parallel_mod

  use mpi
  use flogger
  use const_mod
  use latlon_mesh_mod
  use latlon_halo_mod
  use latlon_parallel_types_mod
  use latlon_parallel_zonal_mod
  use latlon_field_types_mod

  implicit none

  private

  public proc
  public fill_halo
  public zonal_sum
  public zonal_max
  public zonal_avg
  public global_sum
  public global_max

  interface fill_halo
    module procedure fill_halo_2d
    module procedure fill_halo_3d
    module procedure fill_halo_4d
  end interface fill_halo

  interface global_sum
    module procedure global_sum_0d_r4
    module procedure global_sum_0d_r8
  end interface global_sum

  interface global_max
    module procedure global_max_0d_r4
    module procedure global_max_0d_r8
    module procedure global_max_0d_i4
  end interface global_max

contains

  subroutine fill_halo_2d(field, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field2d_type), intent(in) :: field
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, i, j, js, je, nx, mx, hx, hy, ierr
    integer send_req, recv_req
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw)

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)
    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%full_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (west_halo_opt) then
      call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_2d(t1,t2), field%halo(east)%proc_id, 21, &
                        field%d, 1, field%halo(west)%recv_type_2d(t1,t2), field%halo(west)%proc_id, 21, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (east_halo_opt) then
      call MPI_SENDRECV(field%d, 1, field%halo(west)%send_type_2d(t1,t2), field%halo(west)%proc_id, 22, &
                        field%d, 1, field%halo(east)%recv_type_2d(t1,t2), field%halo(east)%proc_id, 22, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_2d(t1,t2), field%halo(north)%proc_id, 23, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_2d(t1,t2), field%halo(south)%proc_id, 23, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (north_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_2d(t1,t2), field%halo(south)%proc_id, 24, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_2d(t1,t2), field%halo(north)%proc_id, 24, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
      call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_2d(t1,t2), field%halo(south)%proc_id, 25, &
                        field%d, 1, field%halo(south)%recv_type_2d(t1,t2), field%halo(south)%proc_id, 25, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,js:0)
      if (field%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = js, 0
          field%d(   1:mx,j) = tmp(hx+1+mx:hx+nx,hy+js-j)
          field%d(mx+1:nx,j) = tmp(hx+1   :hx+mx,hy+js-j)
        end do
      else
        do j = js, 0
          field%d(:,j) = tmp(:,hy+js-j)
        end do
      end if
    end if

    if (north_halo_opt .and. proc%at_north_pole .and. field%halo_cross_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_2d(t1,t2), field%halo(north)%proc_id, 26, &
                        field%d, 1, field%halo(north)%recv_type_2d(t1,t2), field%halo(north)%proc_id, 26, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,je-hy+1:je)
      if (field%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = je - hy + 1, je
          field%d(   1:mx,j) = tmp(hx+1+mx:hx+nx,je+1-j)
          field%d(mx+1:nx,j) = tmp(hx+1   :hx+mx,je+1-j)
        end do
      else
        do j = je - hy + 1, je
          field%d(:,j) = tmp(:,je+1-j)
        end do
      end if
    end if

  end subroutine fill_halo_2d

  subroutine fill_halo_3d(field, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field3d_type), intent(in) :: field
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, t3, i, j, js, je, nx, mx, hx, hy, ierr
    integer send_req, recv_req
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)
    t3 = merge(1, 2, field%full_lev)
    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%half_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (west_halo_opt) then
      call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
                        field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (east_halo_opt) then
      call MPI_SENDRECV(field%d, 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 32, &
                        field%d, 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 32, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (north_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
      call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
                        field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,js:0,:)
      if (field%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = js, 0
          field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
          field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,hy+js-j,:)
        end do
      else
        do j = js, 0
          field%d(:,j,:) = tmp(:,hy+js-j,:)
        end do
      end if
    end if

    if (north_halo_opt .and. proc%at_north_pole .and. field%halo_cross_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
                        field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,je-hy+1:je,:)
      if (field%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = je - hy + 1, je
          field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,je+1-j,:)
          field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,je+1-j,:)
        end do
      else
        do j = je - hy + 1, je
          field%d(:,j,:) = tmp(:,je+1-j,:)
        end do
      end if
    end if

  end subroutine fill_halo_3d

  subroutine fill_halo_4d(field, i4, west_halo, east_halo, south_halo, north_halo)

    type(latlon_field4d_type), intent(in) :: field
    integer, intent(in) :: i4
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
    integer t1, t2, t3, i, j, js, je, nx, mx, hx, hy, ierr
    integer send_req, recv_req
    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))

    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo

    t1 = merge(1, 2, field%full_lon)
    t2 = merge(1, 2, field%full_lat)
    t3 = merge(1, 2, field%full_lev)
    hx = field%halo(1)%lon_hw
    hy = field%halo(1)%lat_hw
    if (field%full_lon) then
      nx = field%mesh%full_nlon
      mx = field%mesh%full_nlon / 2
    else
      nx = field%mesh%half_nlon
      mx = field%mesh%half_nlon / 2
    end if
    if (field%full_lat) then
      js = field%mesh%full_jms
      je = field%mesh%full_jme
    else
      js = field%mesh%half_jms
      je = field%mesh%half_jme
    end if

    if (west_halo_opt) then
      call MPI_SENDRECV(field%d(:,:,:,i4), 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 41, &
                        field%d(:,:,:,i4), 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 41, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (east_halo_opt) then
      call MPI_SENDRECV(field%d(:,:,:,i4), 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 42, &
                        field%d(:,:,:,i4), 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 42, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 43, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 43, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (north_halo_opt) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(field%d(:,:,:,i4), 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 44, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(field%d(:,:,:,i4), 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 44, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
      call MPI_SENDRECV(field%d(:,:,:,i4), 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 45, &
                        field%d(:,:,:,i4), 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 45, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,js:0,:,i4)
      if (field%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = js, 0
          field%d(   1:mx,j,:,i4) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
          field%d(mx+1:nx,j,:,i4) = tmp(hx+1   :hx+mx,hy+js-j,:)
        end do
      else
        do j = js, 0
          field%d(:,j,:,i4) = tmp(:,hy+js-j,:)
        end do
      end if
    end if

    if (north_halo_opt .and. proc%at_north_pole .and. field%halo_cross_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(field%d(:,:,:,i4), 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 46, &
                        field%d(:,:,:,i4), 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 46, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp = field%d(:,je-hy+1:je,:,i4)
      if (field%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = je - hy + 1, je
          field%d(   1:mx,j,:,i4) = tmp(hx+1+mx:hx+nx,je+1-j,:)
          field%d(mx+1:nx,j,:,i4) = tmp(hx+1   :hx+mx,je+1-j,:)
        end do
      else
        do j = je - hy + 1, je
          field%d(:,j,:,i4) = tmp(:,je+1-j,:)
        end do
      end if
    end if

  end subroutine fill_halo_4d

  subroutine global_sum_0d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r4

  subroutine global_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r8

  subroutine global_max_0d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_r4

  subroutine global_max_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_r8

  subroutine global_max_0d_i4(comm, value)

    integer, intent(in) :: comm
    integer(4), intent(inout) :: value

    integer ierr
    integer(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_INT, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_i4

end module latlon_parallel_mod
