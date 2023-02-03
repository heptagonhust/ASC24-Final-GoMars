module parallel_mod

  use mpi
  use flogger
  use const_mod
  use mesh_mod
  use halo_mod
  use process_mod
  use parallel_zonal_mod

  implicit none

  private

  public proc
  public process_init
  public process_barrier
  public process_stop
  public process_final
  public fill_halo
  public zonal_sum
  public zonal_max
  public global_sum
  public global_max

  interface fill_halo
    module procedure fill_halo_2d_r4
    module procedure fill_halo_2d_r8
    module procedure fill_halo_3d_r4
    module procedure fill_halo_3d_r8
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

  subroutine fill_halo_2d_r4(halo, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(halo_type), intent(in) :: halo(:)
    real(4), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(4) tmp(size(array,1),halo(1)%lat_hw)

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = halo(1)%lon_hw
    hy = halo(1)%lat_hw

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, halo(east)%send_type_2d(t1,t2), halo(east)%proc_id, 3, &
                        array, 1, halo(west)%recv_type_2d(t1,t2), halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, halo(west)%send_type_2d(t1,t2), halo(west)%proc_id, 7, &
                        array, 1, halo(east)%recv_type_2d(t1,t2), halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, halo(north)%send_type_2d(t1,t2), halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, halo(south)%recv_type_2d(t1,t2), halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, halo(south)%send_type_2d(t1,t2), halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, halo(north)%recv_type_2d(t1,t2), halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, halo(south)%send_type_2d(t1,t2), halo(south)%proc_id, 11, &
                        array, 1, halo(south)%recv_type_2d(t1,t2), halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,1:hy)
      if (halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1) = tmp(mx+1:nx-hx,j)
          array(mx+1:nx-hx,hy-j+1) = tmp(hx+1:mx   ,j)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1) = tmp(:,j)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, halo(north)%send_type_2d(t1,t2), halo(north)%proc_id, 12, &
                        array, 1, halo(north)%recv_type_2d(t1,t2), halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,ny-hy+1:ny)
      if (halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j) = tmp(mx+1:nx-hx,hy-j+1)
          array(mx+1:nx-hx,ny-hy+j) = tmp(hx+1:mx   ,hy-j+1)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j) = tmp(:,hy-j+1)
        end do
      end if
    end if

  end subroutine fill_halo_2d_r4

  subroutine fill_halo_2d_r8(halo, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(halo_type), intent(in) :: halo(:)
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(8) tmp(size(array,1),halo(1)%lat_hw)

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = halo(1)%lon_hw
    hy = halo(1)%lat_hw

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, halo(east)%send_type_2d(t1,t2), halo(east)%proc_id, 3, &
                        array, 1, halo(west)%recv_type_2d(t1,t2), halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, halo(west)%send_type_2d(t1,t2), halo(west)%proc_id, 7, &
                        array, 1, halo(east)%recv_type_2d(t1,t2), halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, halo(north)%send_type_2d(t1,t2), halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, halo(south)%recv_type_2d(t1,t2), halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, halo(south)%send_type_2d(t1,t2), halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, halo(north)%recv_type_2d(t1,t2), halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, halo(south)%send_type_2d(t1,t2), halo(south)%proc_id, 11, &
                        array, 1, halo(south)%recv_type_2d(t1,t2), halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,1:hy)
      if (halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1) = tmp(mx+1:nx-hx,j)
          array(mx+1:nx-hx,hy-j+1) = tmp(hx+1:mx   ,j)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1) = tmp(:,j)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, halo(north)%send_type_2d(t1,t2), halo(north)%proc_id, 12, &
                        array, 1, halo(north)%recv_type_2d(t1,t2), halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,ny-hy+1:ny)
      if (halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j) = tmp(mx+1:nx-hx,hy-j+1)
          array(mx+1:nx-hx,ny-hy+j) = tmp(hx+1:mx   ,hy-j+1)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j) = tmp(:,hy-j+1)
        end do
      end if
    end if

  end subroutine fill_halo_2d_r8

  subroutine fill_halo_3d_r4(halo, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(halo_type), intent(in) :: halo(:)
    real(4), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, t3, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(4) tmp(size(array,1),halo(1)%lat_hw,size(array,3))

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    t3 = merge(1, 2, merge(full_lev, .true., present(full_lev)))
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = halo(1)%lon_hw
    hy = halo(1)%lat_hw

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, halo(east)%send_type_3d(t1,t2,t3), halo(east)%proc_id, 3, &
                        array, 1, halo(west)%recv_type_3d(t1,t2,t3), halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, halo(west)%send_type_3d(t1,t2,t3), halo(west)%proc_id, 7, &
                        array, 1, halo(east)%recv_type_3d(t1,t2,t3), halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, halo(north)%send_type_3d(t1,t2,t3), halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, halo(south)%recv_type_3d(t1,t2,t3), halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, halo(south)%send_type_3d(t1,t2,t3), halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, halo(north)%recv_type_3d(t1,t2,t3), halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, halo(south)%send_type_3d(t1,t2,t3), halo(south)%proc_id, 11, &
                        array, 1, halo(south)%recv_type_3d(t1,t2,t3), halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:,:) = array(:,1:hy,:)
      if (halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1,:) = tmp(mx+1:nx-hx,j,:)
          array(mx+1:nx-hx,hy-j+1,:) = tmp(hx+1:mx   ,j,:)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1,:) = tmp(:,j,:)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, halo(north)%send_type_3d(t1,t2,t3), halo(north)%proc_id, 12, &
                        array, 1, halo(north)%recv_type_3d(t1,t2,t3), halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      tmp(:,:,:) = array(:,ny-hy+1:ny,:)
      if (halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j,:) = tmp(mx+1:nx-hx,hy-j+1,:)
          array(mx+1:nx-hx,ny-hy+j,:) = tmp(hx+1:mx   ,hy-j+1,:)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j,:) = tmp(:,hy-j+1,:)
        end do
      end if
    end if

  end subroutine fill_halo_3d_r4

  subroutine fill_halo_3d_r8(halo, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(halo_type), intent(in) :: halo(:)
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, t3, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(8) tmp(size(array,1),halo(1)%lat_hw,size(array,3))

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    t3 = merge(1, 2, merge(full_lev, .true., present(full_lev)))
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = halo(1)%lon_hw
    hy = halo(1)%lat_hw

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, halo(east)%send_type_3d(t1,t2,t3), halo(east)%proc_id, 3, &
                        array, 1, halo(west)%recv_type_3d(t1,t2,t3), halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, halo(west)%send_type_3d(t1,t2,t3), halo(west)%proc_id, 7, &
                        array, 1, halo(east)%recv_type_3d(t1,t2,t3), halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, halo(north)%send_type_3d(t1,t2,t3), halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, halo(south)%recv_type_3d(t1,t2,t3), halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, halo(south)%send_type_3d(t1,t2,t3), halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, halo(north)%recv_type_3d(t1,t2,t3), halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, halo(south)%send_type_3d(t1,t2,t3), halo(south)%proc_id, 11, &
                        array, 1, halo(south)%recv_type_3d(t1,t2,t3), halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:,:) = array(:,1:hy,:)
      if (halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1,:) = tmp(mx+1:nx-hx,j,:)
          array(mx+1:nx-hx,hy-j+1,:) = tmp(hx+1:mx   ,j,:)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1,:) = tmp(:,j,:)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, halo(north)%send_type_3d(t1,t2,t3), halo(north)%proc_id, 12, &
                        array, 1, halo(north)%recv_type_3d(t1,t2,t3), halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      tmp(:,:,:) = array(:,ny-hy+1:ny,:)
      if (halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j,:) = tmp(mx+1:nx-hx,hy-j+1,:)
          array(mx+1:nx-hx,ny-hy+j,:) = tmp(hx+1:mx   ,hy-j+1,:)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j,:) = tmp(:,hy-j+1,:)
        end do
      end if
    end if

  end subroutine fill_halo_3d_r8

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

end module parallel_mod
