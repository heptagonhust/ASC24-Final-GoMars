! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_objects_mod

  use simple_physics_types_mod

  implicit none

  type simple_physics_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(simple_state_type) state
    type(simple_tend_type ) tend
  end type simple_physics_objects_type

  type(simple_physics_objects_type), allocatable :: objects(:)

contains

  subroutine simple_physics_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call simple_physics_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state%init(objects(iblk)%mesh)
      call objects(iblk)%tend %init(objects(iblk)%mesh)
    end do

  end subroutine simple_physics_objects_init

  subroutine simple_physics_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state%clear()
        call objects(iblk)%tend %clear()
      end do
    end if

  end subroutine simple_physics_objects_final

end module simple_physics_objects_mod
