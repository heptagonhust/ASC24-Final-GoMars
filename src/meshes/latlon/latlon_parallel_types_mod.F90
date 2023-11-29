! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_parallel_types_mod

  use mpi
  use math_mod, only: round_robin
  use latlon_mesh_mod

  implicit none

  private

  public zonal_circle_type
  public process_type
  public proc

  integer, public :: cart_dim_lon = 1
  integer, public :: cart_dim_lat = 2
  integer, public, parameter :: west  = 1
  integer, public, parameter :: east  = 2
  integer, public, parameter :: south = 3
  integer, public, parameter :: north = 4

  type zonal_circle_type
    integer :: group       = MPI_GROUP_NULL
    integer :: comm        = MPI_COMM_NULL
    integer :: np          = 0
    integer :: id          = MPI_PROC_NULL
    integer :: west_ngb_id = MPI_PROC_NULL
    integer :: east_ngb_id = MPI_PROC_NULL
    integer, allocatable :: recv_type_i4(:,:) ! 0: one level, 1: full_lev, 2: half_lev
    integer, allocatable :: recv_type_r4(:,:) ! 0: one level, 1: full_lev, 2: half_lev
    integer, allocatable :: recv_type_r8(:,:) ! 0: one level, 1: full_lev, 2: half_lev
  contains
    procedure :: init  => zonal_circle_init
    procedure :: clear => zonal_circle_clear
    final :: zonal_circle_final
  end type zonal_circle_type

  type process_neighbor_type
    integer :: id      = MPI_PROC_NULL
    integer :: cart_id = MPI_PROC_NULL
    integer :: orient  = 0
    integer :: lon_hw  = 0
  end type process_neighbor_type

  type process_type
    integer :: comm           = MPI_COMM_NULL
    integer :: cart_comm      = MPI_COMM_NULL
    integer :: group          = MPI_GROUP_NULL
    integer :: cart_group     = MPI_GROUP_NULL
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    integer :: id             = MPI_PROC_NULL          ! MPI process ID
    integer :: cart_id        = MPI_PROC_NULL          ! MPI process ID in cart_comm
    integer idom                                       ! Nest domain index (root domain is 1)
    integer np
    integer nlon
    integer nlat
    integer ids
    integer ide
    integer jds
    integer jde
    logical :: at_south_pole = .false.
    logical :: at_north_pole = .false.
    type(zonal_circle_type) zonal_circle
    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes
    ! Decomposition information
    integer, allocatable :: grid_proc_idmap(:,:) ! Map of grid points to process IDs
    integer, allocatable :: global_grid_id(:,:)  ! Global grid point IDs
    integer, allocatable :: local_grid_id(:,:)   ! Local grid point IDs
    character(30) :: hostname = ''

    integer decomp_type
    integer decomp_loc
  contains
    procedure :: is_root => process_is_root
  end type process_type

  type(process_type) proc

contains

  subroutine zonal_circle_init(this, proc)

    class(zonal_circle_type), intent(inout) :: this
    type(process_type), intent(in) :: proc

    integer ierr, i, nlon, ibeg, iend
    integer cart_coords(2), west_cart_id, east_cart_id, tmp_id(1)
    integer, allocatable :: zonal_proc_id(:)

    allocate(zonal_proc_id(proc%cart_dims(cart_dim_lon)))
    do i = 1, proc%cart_dims(cart_dim_lon)
      cart_coords = proc%cart_coords
      cart_coords(cart_dim_lon) = i - 1
      call MPI_CART_RANK(proc%cart_comm, cart_coords, zonal_proc_id(i), ierr)
    end do
    call MPI_GROUP_INCL(proc%cart_group, size(zonal_proc_id), zonal_proc_id, this%group, ierr)
    call MPI_COMM_CREATE_GROUP(proc%cart_comm, this%group, sum(zonal_proc_id), this%comm, ierr)
    call MPI_COMM_SIZE(this%comm, this%np, ierr)
    call MPI_COMM_RANK(this%comm, this%id, ierr)
    deallocate(zonal_proc_id)

    ! Get IDs of the west and east neighbors in zonal circle comm.
    call MPI_CART_SHIFT(proc%cart_comm, cart_dim_lon-1, 1, west_cart_id, east_cart_id, ierr)
    call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [west_cart_id], this%group, tmp_id, ierr); this%west_ngb_id = tmp_id(1)
    call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [east_cart_id], this%group, tmp_id, ierr); this%east_ngb_id = tmp_id(1)

    if (this%id == 0) then
      ! Integer
      allocate(this%recv_type_i4(this%np,0:2))
      do i = 1, this%np
        nlon = global_mesh%full_nlon
        call round_robin(this%np, i - 1, nlon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(1, [global_mesh%full_nlon], &
                                         [                 nlon], &
                                         [ibeg-1], MPI_ORDER_FORTRAN, MPI_INT, &
                                         this%recv_type_i4(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_i4(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%full_nlev], &
                                         [                 nlon,global_mesh%full_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_INT, &
                                         this%recv_type_i4(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_i4(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%half_nlev], &
                                         [                 nlon,global_mesh%half_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_INT, &
                                         this%recv_type_i4(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_i4(i,2), ierr)
      end do
      ! Single precision
      allocate(this%recv_type_r4(this%np,0:2))
      do i = 1, this%np
        nlon = global_mesh%full_nlon
        call round_robin(this%np, i - 1, nlon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(1, [global_mesh%full_nlon], &
                                         [                 nlon], &
                                         [ibeg-1], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%full_nlev], &
                                         [                 nlon,global_mesh%full_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%half_nlev], &
                                         [                 nlon,global_mesh%half_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,2), ierr)
      end do
      ! Double precision
      allocate(this%recv_type_r8(this%np,0:2))
      do i = 1, this%np
        nlon = global_mesh%full_nlon
        call round_robin(this%np, i - 1, nlon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(1, [global_mesh%full_nlon], &
                                         [                 nlon], &
                                         [ibeg-1], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%full_nlev], &
                                         [                 nlon,global_mesh%full_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%full_nlon,global_mesh%half_nlev], &
                                         [                 nlon,global_mesh%half_nlev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,2), ierr)
      end do
    end if

  end subroutine zonal_circle_init

  subroutine zonal_circle_clear(this)

    class(zonal_circle_type), intent(inout) :: this

    integer i, k, ierr

    if (allocated(this%recv_type_i4)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_i4(i,k), ierr)
        end do
        deallocate(this%recv_type_i4)
      end do
    end if

    if (allocated(this%recv_type_r4)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r4(i,k), ierr)
        end do
        deallocate(this%recv_type_r4)
      end do
    end if

    if (allocated(this%recv_type_r8)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r8(i,k), ierr)
        end do
        deallocate(this%recv_type_r8)
      end do
    end if

    if (this%group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(this%group, ierr)

  end subroutine zonal_circle_clear

  subroutine zonal_circle_final(this)

    type(zonal_circle_type), intent(inout) :: this

    call this%clear()

  end subroutine zonal_circle_final

  pure logical function process_is_root(this) result(res)

    class(process_type), intent(in) :: this

    res = this%id == 0

  end function process_is_root

end module latlon_parallel_types_mod
