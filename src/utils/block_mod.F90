module block_mod

  use mpi
  use flogger
  use namelist_mod
  use mesh_mod
  use dynamics_types_mod
  use physics_types_mod
  use adv_batch_mod
  use filter_types_mod
  use halo_mod
  use allocator_mod

  implicit none

  private

  public block_type
  public blocks
  public global_mesh
  public mesh_type
  public static_type
  public dstate_type
  public pstate_type
  public dtend_type
  public ptend_type

  type block_type
    integer id
    type(mesh_type) filter_mesh
    type(mesh_type) mesh
    type(static_type) static
    type(dstate_type), allocatable :: dstate(:)
    type(pstate_type) pstate
    type(dtend_type), allocatable :: dtend(:)
    type(ptend_type) ptend
    type(adv_batch_type) adv_batch_pt
    type(adv_batch_type), allocatable :: adv_batches(:)
    type(filter_type) big_filter
    type(filter_type) small_filter
    type(halo_type), allocatable :: filter_halo(:)
    type(halo_type), allocatable :: halo(:)
  contains
    procedure :: init_stage_1 => block_init_stage_1
    procedure :: init_stage_2 => block_init_stage_2
    procedure :: clear => block_clear
    final :: block_final
  end type block_type

  type(block_type), allocatable :: blocks(:)

contains

  subroutine block_init_stage_1(this, id, ids, ide, jds, jde)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: ids
    integer, intent(in) :: ide
    integer, intent(in) :: jds
    integer, intent(in) :: jde

    this%id = id

    call this%filter_mesh%init_from_parent(global_mesh, this%id, ids, ide, jds, jde)
    call this%mesh%init_from_parent(global_mesh, this%id, ids, ide, jds, jde)
    call this%big_filter%init(this%filter_mesh, 'big_filter')
    call this%small_filter%init(this%filter_mesh, 'small_filter')

  end subroutine block_init_stage_1

  subroutine block_init_stage_2(this)

    class(block_type), intent(inout) :: this

    integer i

    call this%filter_mesh%reinit()

    if (.not. allocated(this%dstate)) then
      select case (trim(time_scheme))
      case ('euler')
        allocate(this%dstate(2))
        allocate(this%dtend (2))
      case ('pc2', 'wrfrk3')
        allocate(this%dstate(3))
        allocate(this%dtend (3))
      case default
        if (this%id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      do i = 1, size(this%dstate)
        call this%dstate(i)%init(this%filter_mesh, this%mesh)
      end do
      do i = 1, size(this%dtend)
        call this%dtend(i)%init(this%filter_mesh, this%mesh)
      end do
      call this%static%init(this%filter_mesh, this%mesh)
      call this%pstate%init(this%mesh)
      call this%ptend%init(this%mesh)
    end if

  end subroutine block_init_stage_2

  subroutine block_clear(this)

    class(block_type), intent(inout) :: this

    integer i

    call this%filter_mesh%clear()
    call this%mesh%clear()
    call this%big_filter%clear()
    call this%small_filter%clear()
    do i = 1, size(this%dstate)
      call this%dstate(i)%clear()
    end do
    do i = 1, size(this%dtend)
      call this%dtend(i)%clear()
    end do
    call this%pstate%clear()
    call this%ptend %clear()
    call this%adv_batch_pt%clear()
    if (allocated(this%adv_batches)) then
      do i = 1, size(this%adv_batches)
        call this%adv_batches(i)%clear()
      end do
    end if
    if (allocated(this%halo)) then
      do i = 1, size(this%halo)
        call this%halo(i)%clear()
      end do
    end if

    if (allocated(this%dstate)) deallocate(this%dstate)
    if (allocated(this%dtend)) deallocate(this%dtend)
    if (allocated(this%adv_batches)) deallocate(this%adv_batches)
    if (allocated(this%halo)) deallocate(this%halo)

  end subroutine block_clear

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    call this%clear()

  end subroutine block_final

end module block_mod
