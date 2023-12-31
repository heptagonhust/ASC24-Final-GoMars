module wrf_physics_types_mod

  use const_mod
  use physics_types_mod
  use wrf_namelist_mod
  use noahlsm_types_mod

  implicit none

  private

  type, extends(physics_state_type) :: wrf_state_type
    type(noahlsm_state_type) noahlsm
  contains
    procedure :: init  => wrf_state_init
    procedure :: clear => wrf_state_clear
    final wrf_state_final
  end type wrf_state_type

contains

  subroutine wrf_state_init(this, mesh)

    class(wrf_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    call this%noahlsm%init(mesh)

    call this%physics_state_init(mesh)

  end subroutine wrf_state_init

  subroutine wrf_state_clear(this)

    class(wrf_state_type), intent(inout) :: this

    call this%noahlsm%clear()

    call this%physics_state_clear()

  end subroutine wrf_state_clear

  subroutine wrf_state_final(this)

    type(wrf_state_type), intent(inout) :: this

    call this%clear()

  end subroutine wrf_state_final

end module wrf_physics_types_mod
