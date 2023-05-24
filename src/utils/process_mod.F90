module process_mod

  use mpi
  use flogger
  use string
  use const_mod
  use namelist_mod
  use latlon_mesh_mod
  use block_mod
  use math_mod
  use perf_mod
  use parallel_types_mod

  implicit none

  private

  public process_init
  public process_create_blocks
  public process_barrier
  public process_stop
  public process_final
  public proc
  public zonal_circle_type

  integer, public, parameter :: decomp_2d_simple     = 1

  integer, public, parameter :: decomp_normal_region = 5

contains

  subroutine process_init(comm)

    integer, intent(in), optional :: comm

    integer ierr, n

    if (present(comm)) then
      proc%comm = comm
    else
      call MPI_INIT(ierr)
      proc%comm = MPI_COMM_WORLD
    end if
    call perf_init()
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)
    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)
    call MPI_GET_PROCESSOR_NAME(proc%hostname, n, ierr)

    call setup_mpi_simple()
    call decompose_domains()
    call setup_zonal_comm()

    select case (proc_layout)
    case ('lon>lat')
      cart_dim_lon = 2
      cart_dim_lat = 1
      call log_notice('Process layout: lon > lat', pid=proc%id)
    case ('lat>lon')
      cart_dim_lon = 1
      cart_dim_lat = 2
      call log_notice('Process layout: lat > lon', pid=proc%id)
    case default
      call log_error('Invalid proc_layout ' // trim(proc_layout) // '!', pid=proc%id)
    end select

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

  subroutine setup_mpi_simple()

    integer ierr, np, tmp_comm, i
    logical periods(2)

    proc%decomp_type = decomp_2d_simple
    proc%decomp_loc  = decomp_normal_region

    if (nproc_lon(1) * nproc_lat(1) /= proc%np) then
      nproc_lat(1) = global_mesh%full_nlat / 2
      nproc_lon(1) = proc%np / nproc_lat(1)
    end if
    call log_notice('Set process layout to ' // to_str(nproc_lon(1)) // ' x ' // to_str(nproc_lat(1)) // '.', pid=proc%id)

    if (nproc_lon(1) * nproc_lat(1) == proc%np) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np = 0
      do i = 1, 1
        np = np + nproc_lon(i) * nproc_lat(i)
      end do
      if (proc%np /= np .and. proc%is_root()) then
        call log_notice('Namelist nproc_lon and nproc_lat are not compatible with MPI runtime. Reset to MPI runtime.')
        nproc_lat(1) = proc%np
      end if
      ! Set the process topology into proc object.
      np = 0
      do i = 1, 1
        np = np + nproc_lon(i) * nproc_lat(i)
        if (proc%id + 1 <= np) then
          proc%cart_dims(cart_dim_lon) = nproc_lon(i)
          proc%cart_dims(cart_dim_lat) = nproc_lat(i)
          proc%idom = i
          exit
        end if
      end do
    else
      proc%cart_dims = [merge(1, proc%np, cart_dim_lon == 1), merge(proc%np, 1, cart_dim_lat == 2)]
      proc%idom = 1
    end if
    ! Check decomposition dimensions.
    if (proc%cart_dims(cart_dim_lon) /= 1 .and. mod(proc%cart_dims(cart_dim_lon), 2) /= 0) then
      call process_stop(1, 'nproc_lon should be an even number!')
    end if
    ! Set MPI process topology.
    periods = [cart_dim_lon==1,cart_dim_lon==2]
    call MPI_COMM_SPLIT(proc%comm, proc%idom, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%cart_comm, ierr)
    call MPI_COMM_GROUP(proc%cart_comm, proc%cart_group, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_RANK(proc%cart_comm, proc%cart_id, ierr)
    call MPI_CART_COORDS(proc%cart_comm, proc%cart_id, 2, proc%cart_coords, ierr)

  end subroutine setup_mpi_simple

  subroutine decompose_domains()

    integer ierr, tmp_id(1), i, j, ip, global_ig, local_ig
    integer, allocatable :: all_ids(:), all_ide(:), all_jds(:), all_jde(:)
    logical correct, all_correct

    ! Set neighborhood of the process.
    if (allocated(proc%ngb)) deallocate(proc%ngb)
    select case (proc%decomp_loc)
    case (decomp_normal_region)
      allocate(proc%ngb(4))
      call MPI_CART_SHIFT(proc%cart_comm, cart_dim_lon-1, 1, proc%ngb(west )%cart_id, proc%ngb(east )%cart_id, ierr)
      call MPI_CART_SHIFT(proc%cart_comm, cart_dim_lat-1, 1, proc%ngb(south)%cart_id, proc%ngb(north)%cart_id, ierr)
    end select

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, size(proc%ngb)
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    ! Handle processes at poles.
    if (proc%ngb(south)%id == MPI_PROC_NULL) then
      if (cart_dim_lon == 1) then
        i = proc%id + proc%cart_dims(cart_dim_lon) / 2 * proc%cart_dims(cart_dim_lat)
        if (i >= proc%np) i = i - proc%np
      else
        i = proc%id + proc%cart_dims(cart_dim_lon) / 2
        if (i >= proc%cart_dims(cart_dim_lon)) i = i - proc%cart_dims(cart_dim_lon)
      end if
      proc%ngb(south)%id = i
      proc%at_south_pole = .true.
    end if
    if (proc%ngb(north)%id == MPI_PROC_NULL) then
      if (cart_dim_lon == 1) then
        i = proc%id + proc%cart_dims(cart_dim_lon) / 2 * proc%cart_dims(cart_dim_lat)
        if (i >= proc%np) i = i - proc%np
      else
        i = proc%id + proc%cart_dims(cart_dim_lat) / 2
        if (i >= proc%np) i = i - proc%cart_dims(cart_dim_lon)
      end if
      proc%ngb(north)%id = i
      proc%at_north_pole = .true.
    end if

    ! Set initial values for nlon, nlat, ids, jds.
    proc%nlon = global_mesh%full_nlon
    select case (proc%decomp_loc)
    case (decomp_normal_region)
      proc%nlat = global_mesh%full_nlat
    end select

    call round_robin(proc%cart_dims(cart_dim_lon), proc%cart_coords(cart_dim_lon), proc%nlon, proc%ids, proc%ide)
    call round_robin(proc%cart_dims(cart_dim_lat), proc%cart_coords(cart_dim_lat), proc%nlat, proc%jds, proc%jde)

    correct = .true.
    if (proc%nlat < 2) then
      correct = .false.
    end if
    call MPI_ALLREDUCE(correct, all_correct, 1, MPI_LOGICAL, MPI_LAND, proc%comm, ierr)
    if (.not. all_correct) then
      call log_error('Decomposed grid size along latitude should be >= 2!', __FILE__, __LINE__, proc%id)
    end if

    ! Set grid_proc_idmap for later use.
    allocate(proc%grid_proc_idmap(nlon,nlat))
    allocate(proc%global_grid_id (nlon,nlat))
    allocate(proc%local_grid_id  (nlon,nlat))
    allocate(all_ids(proc%np), all_ide(proc%np))
    allocate(all_jds(proc%np), all_jde(proc%np))
    call MPI_ALLGATHER(proc%ids, 1, MPI_INTEGER, all_ids, 1, MPI_INTEGER, proc%cart_comm, ierr)
    call MPI_ALLGATHER(proc%ide, 1, MPI_INTEGER, all_ide, 1, MPI_INTEGER, proc%cart_comm, ierr)
    call MPI_ALLGATHER(proc%jds, 1, MPI_INTEGER, all_jds, 1, MPI_INTEGER, proc%cart_comm, ierr)
    call MPI_ALLGATHER(proc%jde, 1, MPI_INTEGER, all_jde, 1, MPI_INTEGER, proc%cart_comm, ierr)
    global_ig = 1
    do ip = 1, proc%np
      proc%grid_proc_idmap(all_ids(ip):all_ide(ip), all_jds(ip):all_jde(ip)) = ip
      local_ig = 1
      do j = all_jds(ip), all_jde(ip)
        do i = all_ids(ip), all_ide(ip)
          proc%global_grid_id(i,j) = global_ig
          global_ig = global_ig + 1
          proc%local_grid_id(i,j) = local_ig
          local_ig = local_ig + 1
        end do
      end do
    end do
    deallocate(all_ids, all_ide, all_jds, all_jde)

  end subroutine decompose_domains

  subroutine setup_zonal_comm()

    if (proc%idom == 1) then ! Only root domain has polar region.
      call proc%zonal_circle%init(proc)
    end if

  end subroutine setup_zonal_comm

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
