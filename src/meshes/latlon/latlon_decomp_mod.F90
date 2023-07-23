module latlon_decomp_mod

  use mpi
  use string
  use flogger
  use math_mod, only: round_robin
  use latlon_mesh_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public latlon_decomp_run

  integer, public, parameter :: decomp_2d_simple     = 1
  integer, public, parameter :: decomp_normal_region = 5

contains

  subroutine latlon_decomp_run(proc_layout, nproc_lon, nproc_lat, ierr)

    character(*), intent(in) :: proc_layout
    integer, intent(inout) :: nproc_lon(:)
    integer, intent(inout) :: nproc_lat(:)
    integer, intent(out) :: ierr

    integer np, tmp_comm, tmp_id(1), i, j, n, ip, global_ig, local_ig
    integer, allocatable :: all_ids(:), all_ide(:), all_jds(:), all_jde(:)
    logical periods(2), correct, all_correct

    select case (proc_layout)
    case ('lon>lat')
      cart_dim_lon = 2
      cart_dim_lat = 1
    case ('lat>lon')
      cart_dim_lon = 1
      cart_dim_lat = 2
    end select

    proc%decomp_type = decomp_2d_simple
    proc%decomp_loc  = decomp_normal_region

    if (nproc_lon(1) * nproc_lat(1) /= proc%np) then
      ! User does not set process dimensions in the namelist, so we set them here.
      if (proc%np < global_mesh%full_nlon / 2) then
        nproc_lat(1) = proc%np
      else
        nproc_lat(1) = global_mesh%full_nlat / 3
      end if
      if (mod(proc%np, nproc_lat(1)) /= 0) then
        ierr = 3
        return
      end if
      nproc_lon(1) = proc%np / nproc_lat(1)
    end if

    if (nproc_lon(1) * nproc_lat(1) == proc%np) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np = 0
      do i = 1, 1
        np = np + nproc_lon(i) * nproc_lat(i)
      end do
      if (proc%np /= np .and. proc%is_root()) then
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
    if (proc%is_root()) then
      call log_notice('Process layout is ' // to_str(nproc_lon(1)) // 'x' // to_str(nproc_lat(1)) // '.')
    end if
    ! Check decomposition dimensions.
    if (proc%cart_dims(cart_dim_lon) /= 1 .and. mod(proc%cart_dims(cart_dim_lon), 2) /= 0) then
      ierr = 1 ! nproc_lon should be an even number!
      return
    end if
    ! Set MPI process topology.
    periods = [cart_dim_lon==1,cart_dim_lon==2]
    call MPI_COMM_SPLIT(proc%comm, proc%idom, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%cart_comm, ierr)
    call MPI_COMM_GROUP(proc%cart_comm, proc%cart_group, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_RANK(proc%cart_comm, proc%cart_id, ierr)
    call MPI_CART_COORDS(proc%cart_comm, proc%cart_id, 2, proc%cart_coords, ierr)

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
        i = proc%id + proc%cart_dims(cart_dim_lon) / 2
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
      ierr = 2 ! Decomposed grid size along latitude should be >= 2!
      return
    end if

    ! Set grid_proc_idmap for later use.
    allocate(proc%grid_proc_idmap(global_mesh%full_nlon,global_mesh%full_nlat))
    allocate(proc%global_grid_id (global_mesh%full_nlon,global_mesh%full_nlat))
    allocate(proc%local_grid_id  (global_mesh%full_nlon,global_mesh%full_nlat))
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

    ! Create zonal circle communicator on root domains.
    if (proc%idom == 1) then
      call proc%zonal_circle%init(proc)
    end if

    ierr = 0

  end subroutine latlon_decomp_run

end module latlon_decomp_mod
