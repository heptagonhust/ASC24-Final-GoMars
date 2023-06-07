module process_mod

  use mpi
  use flogger
  use string
  use namelist_mod
  use block_mod
  use math_mod
  use perf_mod
  use latlon_mesh_mod
  use latlon_parallel_types_mod
  use latlon_decomp_mod

  implicit none

  private

  public process_init
  public process_create_blocks
  public process_barrier
  public process_stop
  public process_final
  public proc
  public zonal_circle_type

contains

  subroutine process_init(comm)

    integer, intent(in), optional :: comm

    integer thread, ierr, n

    if (present(comm)) then
      proc%comm = comm
    else
      call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, thread, ierr)
      proc%comm = MPI_COMM_WORLD
      call perf_init()
    end if
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)
    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)
    call MPI_GET_PROCESSOR_NAME(proc%hostname, n, ierr)

    call latlon_decomp_run(proc_layout, nproc_lon, nproc_lat, ierr)

  end subroutine process_init

  subroutine process_barrier()

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)

  end subroutine process_barrier

  subroutine process_stop(code, msg)

    integer, intent(in) :: code
    character(*), intent(in), optional :: msg

    integer ierr

    if (present(msg)) then
      call log_warning(msg, pid=proc%id)
    end if
    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer i, ierr

    do i = 1, size(blocks)
      call blocks(i)%clear()
    end do

    if (allocated(proc%ngb)) deallocate(proc%ngb)
    if (allocated(proc%grid_proc_idmap)) deallocate(proc%grid_proc_idmap)
    if (allocated(proc%global_grid_id)) deallocate(proc%global_grid_id)
    if (allocated(proc%local_grid_id)) deallocate(proc%local_grid_id)
    if (allocated(blocks)) deallocate(blocks)
    if (proc%group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group, ierr)
    if (proc%cart_group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%cart_group, ierr)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

  subroutine process_create_blocks()

    integer lon_hw, i, j, dtype
    integer ierr, status(MPI_STATUS_SIZE)

    if (.not. allocated(blocks)) allocate(blocks(1))

    call blocks(1)%init_stage_1(1, proc%ids, proc%ide, proc%jds, proc%jde)

    ! Each process calculate lon_hw from its big_filter%ngrid_lat(:) and big_filter%ngrid_lon(:).
    lon_hw = global_mesh%lon_hw
    do j = blocks(1)%mesh%half_jds, blocks(1)%mesh%half_jde
      lon_hw = max(lon_hw, (blocks(1)%big_filter%ngrid_lat(j) - 1) / 2)
    end do
    do j = blocks(1)%mesh%full_jds, blocks(1)%mesh%full_jde
      lon_hw = max(lon_hw, (blocks(1)%big_filter%ngrid_lon(j) - 1) / 2)
    end do
    lon_hw = max(lon_hw, global_mesh%lon_hw)

    ! Get lon_hw from southern and northern neighbors.
    if (proc%ngb(south)%id /= MPI_PROC_NULL) then
      call MPI_SENDRECV(lon_hw, 1, MPI_INT, proc%ngb(south)%id, 100, &
                        proc%ngb(south)%lon_hw, 1, MPI_INT, proc%ngb(south)%id, 100, &
                        proc%comm, status, ierr)
    else
      proc%ngb(south)%lon_hw = 0
    end if
    if (proc%ngb(north)%id /= MPI_PROC_NULL) then
      call MPI_SENDRECV(lon_hw, 1, MPI_INT, proc%ngb(north)%id, 100, &
                        proc%ngb(north)%lon_hw, 1, MPI_INT, proc%ngb(north)%id, 100, &
                        proc%comm, status, ierr)
    else
      proc%ngb(north)%lon_hw = 0
    end if

    if (lon_hw > blocks(1)%mesh%full_nlon) then
      call log_error('Too large zonal halo width ' // to_str(lon_hw) // '!', __FILE__, __LINE__)
    end if

    call global_mesh%reinit(lon_hw)
    call blocks(1)%init_stage_2()

    if (proc%is_root()) then
      call log_notice('Maximum zonal halo width is ' // to_str(lon_hw) // '.')
    end if

    select case (r8)
    case (4)
      dtype = MPI_REAL
    case (8)
      dtype = MPI_DOUBLE
    case (16)
      dtype = MPI_REAL16
    case default
      call log_error('Unsupported parameter r8!')
    end select

    ! Setup halos (only normal halos for the time being).
    allocate(blocks(1)%filter_halo(size(proc%ngb)))
    allocate(blocks(1)%halo(size(proc%ngb)))
    do i = 1, size(proc%ngb)
      proc%ngb(i)%orient = i
      select case (proc%ngb(i)%orient)
      case (west, east)
        call blocks(1)%filter_halo(i)%init(blocks(1)%filter_mesh, proc%ngb(i)%orient, dtype,   &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, &
                                    jds=proc%jds, jde=proc%jde)
        call blocks(1)%halo(i)%init(blocks(1)%mesh, proc%ngb(i)%orient, dtype,   &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, &
                                    jds=proc%jds, jde=proc%jde)
      case (south, north)
        lon_hw = 0 ! min(blocks(1)%filter_mesh%lon_hw, proc%ngb(i)%lon_hw)
        call blocks(1)%filter_halo(i)%init(blocks(1)%filter_mesh, proc%ngb(i)%orient, dtype,   &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, &
                                    ids=proc%ids, ide=proc%ide, lon_hw=lon_hw,   &
                                    at_south_pole=proc%at_south_pole,            &
                                    at_north_pole=proc%at_north_pole)
        lon_hw = min(blocks(1)%mesh%lon_hw, proc%ngb(i)%lon_hw)
        call blocks(1)%halo(i)%init(blocks(1)%mesh, proc%ngb(i)%orient, dtype,   &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id, &
                                    ids=proc%ids, ide=proc%ide, lon_hw=lon_hw,   &
                                    at_south_pole=proc%at_south_pole,            &
                                    at_north_pole=proc%at_north_pole)
      end select
    end do

  end subroutine process_create_blocks

end module process_mod
