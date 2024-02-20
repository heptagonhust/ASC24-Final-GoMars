module wrf_objects_mod

  use physics_types_mod
  use wrf_physics_types_mod

  implicit none

  type wrf_objects_type
    type(physics_mesh_type), pointer :: mesh
    type(wrf_state_type) state
    type(wrf_tend_type) tend
  end type wrf_objects_type

  type(wrf_objects_type), allocatable, target :: objects(:)

contains

  subroutine wrf_objects_init(mesh)

    type(physics_mesh_type), intent(in), target :: mesh(:)

    integer nblk, iblk

    call wrf_objects_final()

    nblk = size(mesh)
    allocate(objects(nblk))
    do iblk = 1, nblk
      objects(iblk)%mesh => mesh(iblk)
      call objects(iblk)%state %init(objects(iblk)%mesh)
      call objects(iblk)%tend  %init(objects(iblk)%mesh)
    end do

  end subroutine wrf_objects_init

  subroutine wrf_objects_final()

    integer iblk

    if (allocated(objects)) then
      do iblk = 1, size(objects)
        call objects(iblk)%state %clear()
        call objects(iblk)%tend  %clear()
      end do
    end if

  end subroutine wrf_objects_final

end module wrf_objects_mod
