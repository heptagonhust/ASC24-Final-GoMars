module dyn_comp

  use block_mod

  implicit none

  private

  public dyn_init
  public dyn_final
  public dyn_register
  public dyn_import_t
  public dyn_export_t

  type pstate_ptr_type
    type(pstate_type), pointer :: ptr => null()
  end type pstate_ptr_type

  type dyn_import_t
    type(pstate_ptr_type), allocatable :: pstate(:)
  contains
    final :: dyn_import_final
  end type dyn_import_t

  type ptend_ptr_type
    type(ptend_type), pointer :: ptr => null()
  end type ptend_ptr_type

  type dyn_export_t
    type(ptend_ptr_type), allocatable :: ptend(:)
  contains
    final :: dyn_export_final
  end type dyn_export_t

contains

  subroutine dyn_init(dyn_in, dyn_out)

    type(dyn_import_t), intent(out) :: dyn_in
    type(dyn_export_t), intent(out) :: dyn_out

    integer iblk

    allocate(dyn_in%pstate(size(blocks)))
    allocate(dyn_out%ptend(size(blocks)))
    do iblk = 1, size(blocks)
      dyn_in%pstate(iblk)%ptr => blocks(iblk)%pstate
      dyn_out%ptend(iblk)%ptr => blocks(iblk)%ptend
    end do

  end subroutine dyn_init

  subroutine dyn_final()

  end subroutine dyn_final

  subroutine dyn_register()

  end subroutine dyn_register

  subroutine dyn_import_final(this)

    type(dyn_import_t), intent(inout) :: this

    if (allocated(this%pstate)) deallocate(this%pstate)

  end subroutine dyn_import_final

  subroutine dyn_export_final(this)

    type(dyn_export_t), intent(inout) :: this

    if (allocated(this%ptend)) deallocate(this%ptend)

  end subroutine dyn_export_final

end module dyn_comp