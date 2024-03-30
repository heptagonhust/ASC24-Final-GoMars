! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cubed_sphere_mesh_mod

  use cubed_sphere_panel_mod

  implicit none

  private

  public cubed_sphere_mesh_type
  public cubed_sphere_panel_type

  integer, public, parameter :: grid_loc_cell  = 1
  integer, public, parameter :: grid_loc_vtx   = 2
  integer, public, parameter :: grid_loc_xedge = 3
  integer, public, parameter :: grid_loc_yedge = 4
  integer, public, parameter :: grid_loc_zedge = 5

  type cubed_sphere_mesh_type
    logical :: initialized = .false.
    integer :: grid_type   = 1
    type(cubed_sphere_panel_type), allocatable :: panels(:)
  contains
    procedure :: init      => cubed_sphere_mesh_init
    procedure :: clear     => cubed_sphere_mesh_clear
    final :: cubed_sphere_mesh_final
  end type cubed_sphere_mesh_type

contains

  subroutine cubed_sphere_mesh_init(this, proj_type, nc, nz, hw, ipn, ids, jds, ndx, ndy)

    class(cubed_sphere_mesh_type), intent(inout) :: this
    integer, intent(in) :: proj_type
    integer, intent(in) :: nc             ! Panel grid size along x axis and y axis
    integer, intent(in) :: nz             ! Level size
    integer, intent(in) :: hw             ! Halo width
    integer, intent(in), optional :: ipn  ! Panel index
    integer, intent(in), optional :: ids  ! Domain start index along x axis
    integer, intent(in), optional :: jds  ! Domain start index along y axis
    integer, intent(in), optional :: ndx  ! Domain grid size along x axis
    integer, intent(in), optional :: ndy  ! Domain grid size along y axis

    integer i, il, ir, ib, it

    call this%clear()

    !       Staircase                       Odd Tile                 Even Tile
    !
    !              ___________               _____                     _____
    !             |     |     |             |     |                   |     |
    !             |  5  |  6  |             | n+2 |                   | n+1 |
    !        _____|_____|_____|        _____|_____|_____         _____|_____|_____
    !       |     |     |             |     |     |     |       |     |     |     |
    !       |  3  |  4  |             | n+4 |  n  | n+1 |       | n+5 |  n  | n+2 |
    !  _____|_____|_____|             |_____|_____|_____|       |_____|_____|_____|
    ! |     |     |                         |     |                   |     |
    ! |  1  |  2  |                         | n+5 |                   | n+4 |
    ! |_____|_____|                         |_____|                   |_____|
    !

    allocate(this%panels(6))

    ! Connect panels.
    do i = 1, 6
      if (mod(i, 2) == 0) then
        ! Even tile
        il = i + 5; if (il > 6) il = il - 6
        ir = i + 2; if (ir > 6) ir = ir - 6
        ib = i + 4; if (ib > 6) ib = ib - 6
        it = i + 1; if (it > 6) it = it - 6
      else
        ! Odd tile
        il = i + 4; if (il > 6) il = il - 6
        ir = i + 1; if (ir > 6) ir = ir - 6
        ib = i + 5; if (ib > 6) ib = ib - 6
        it = i + 2; if (it > 6) it = it - 6
      end if
      call this%panels(i)%connect(i, this%panels(il), this%panels(ir), this%panels(ib), this%panels(it))
    end do

    if (present(ipn)) then
      ! Parallel or regional mesh
      call this%panels(ipn)%init(proj_type, nc=nc, nz=nz, hw=hw, ids=ids, jds=jds, ndx=ndx, ndy=ndy, active=.true.)
    else
      ! Global mesh
      do i = 1, 6
        call this%panels(i)%init(proj_type, nc=nc, nz=nz, hw=hw)
        this%panels(i)%active = .true.
      end do
    end if

    this%initialized = .true.

  end subroutine cubed_sphere_mesh_init

  subroutine cubed_sphere_mesh_clear(this)

    class(cubed_sphere_mesh_type), intent(inout) :: this

    if (allocated(this%panels)) deallocate(this%panels)

    this%initialized = .false.

  end subroutine cubed_sphere_mesh_clear

  subroutine cubed_sphere_mesh_final(this)

    type(cubed_sphere_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine cubed_sphere_mesh_final

end module cubed_sphere_mesh_mod
