! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cubed_sphere_array_mod

  use const_mod
  use cubed_sphere_mesh_mod

  implicit none

  private

  public cubed_sphere_array_type

  type cubed_sphere_values_type
    real(r8), allocatable :: values(:,:,:)
  contains
    procedure :: init     => cubed_sphere_values_init
    procedure :: clear    => cubed_sphere_values_clear
    final :: cubed_sphere_values_final
  end type cubed_sphere_values_type

  type cubed_sphere_array_type
    logical :: initialized      = .false.
    character(30) :: name       = ''
    character(30) :: units      = ''
    character(50) :: long_name  = ''
    integer       :: loc        = 0
    type(cubed_sphere_mesh_type), pointer :: mesh => null()
    type(cubed_sphere_values_type), allocatable :: panels(:)
  contains
    procedure :: init     => cubed_sphere_array_init
    procedure :: clear    => cubed_sphere_array_clear
    final :: cubed_sphere_array_final
  end type cubed_sphere_array_type

contains

  subroutine cubed_sphere_values_init(this, panel, loc)

    class(cubed_sphere_values_type), intent(inout) :: this
    type(cubed_sphere_panel_type), intent(in) :: panel
    integer, intent(in) :: loc

    select case (loc)
    case (grid_loc_cell)
      allocate(this%values(panel%full_ims:panel%full_ime,panel%full_jms:panel%full_jme,panel%full_kms:panel%full_kme))
    case (grid_loc_vtx)
      allocate(this%values(panel%half_ims:panel%half_ime,panel%half_jms:panel%half_jme,panel%full_kms:panel%full_kme))
    case (grid_loc_xedge)
      allocate(this%values(panel%half_ims:panel%half_ime,panel%full_jms:panel%full_jme,panel%full_kms:panel%full_kme))
    case (grid_loc_yedge)
      allocate(this%values(panel%full_ims:panel%full_ime,panel%half_jms:panel%half_jme,panel%full_kms:panel%full_kme))
    case (grid_loc_zedge)
      allocate(this%values(panel%full_ims:panel%full_ime,panel%full_jms:panel%full_jme,panel%half_kms:panel%half_kme))
    end select

  end subroutine cubed_sphere_values_init

  subroutine cubed_sphere_values_clear(this)

    class(cubed_sphere_values_type), intent(inout) :: this

    if (allocated(this%values)) deallocate(this%values)

  end subroutine cubed_sphere_values_clear

  subroutine cubed_sphere_values_final(this)

    type(cubed_sphere_values_type), intent(inout) :: this

    call this%clear()

  end subroutine cubed_sphere_values_final

  subroutine cubed_sphere_array_init(this, mesh, name, units, long_name, loc)

    class(cubed_sphere_array_type), intent(inout) :: this
    type(cubed_sphere_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: name
    character(*), intent(in) :: units
    character(*), intent(in) :: long_name
    integer, intent(in) :: loc

    integer i

    call this%clear()

    this%name = name
    this%units = units
    this%long_name = long_name

    allocate(this%panels(size(mesh%panels)))
    do i = 1, size(mesh%panels)
      if (mesh%panels(i)%initialized) then
        call this%panels(i)%init(mesh%panels(i), loc)
      end if
    end do

    this%initialized = .true.

  end subroutine cubed_sphere_array_init

  subroutine cubed_sphere_array_clear(this)

    class(cubed_sphere_array_type), intent(inout) :: this

    integer i

    if (allocated(this%panels)) then
      do i = 1, size(this%panels)
        call this%panels(i)%clear()
      end do
      deallocate(this%panels)
    end if
    this%initialized = .false.

  end subroutine cubed_sphere_array_clear

  subroutine cubed_sphere_array_final(this)

    type(cubed_sphere_array_type), intent(inout) :: this

    call this%clear()

  end subroutine cubed_sphere_array_final

end module cubed_sphere_array_mod