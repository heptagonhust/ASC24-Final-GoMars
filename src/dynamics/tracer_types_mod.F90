module tracer_types_mod

  use const_mod
  use allocator_mod
  use mesh_mod

  implicit none

  integer nbatches
  integer ntracers
  integer ntracers_water

  integer idx_qv
  integer idx_qc
  integer idx_qi
  integer idx_qr
  integer idx_qs
  integer idx_qg
  integer idx_qh
  integer idx_qo3
  integer idx_qso2

  character(32), allocatable :: batch_names(:)
  real(r8), allocatable :: batch_dts(:)
  character(32), allocatable :: tracer_batches(:)
  character(32), allocatable :: tracer_names(:)
  character(32), allocatable :: tracer_long_names(:)
  character(32), allocatable :: tracer_units(:)
  ! Follow other model's definition:
  ! 0 - Generic tracer
  ! 1 - Prognostic chemical tracer
  ! 2 - Diagnostic chemical tracer
  integer, allocatable :: tracer_types(:)

  type tracers_type
    logical :: is_initialized = .false.
    type(mesh_type), pointer :: mesh => null()
    type(mesh_type), pointer :: filter_mesh => null()
    real(r8), allocatable :: q(:,:,:,:)
    ! Some diagnostics:
    real(r8), allocatable :: qm(:,:,:) ! Total moisture or water substances
  contains
    procedure :: init => tracers_init
    procedure :: clear => tracers_clear
    final :: tracers_final
  end type tracers_type

  type(tracers_type), allocatable, target :: tracers(:) ! (blocks)

contains

  subroutine tracers_init(this, mesh, filter_mesh)

    class(tracers_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    type(mesh_type), intent(in), target :: filter_mesh

    call this%clear()

    this%mesh => mesh
    this%filter_mesh => filter_mesh
    call allocate_array(filter_mesh, this%q, extra_dim=ntracers, full_lon=.true., full_lat=.true., full_lev=.true.)
    if (idx_qv > 0) then
      call allocate_array(mesh, this%qm, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

    this%is_initialized = .true.

  end subroutine tracers_init

  subroutine tracers_clear(this)

    class(tracers_type), intent(inout) :: this

    if (allocated(this%q )) deallocate(this%q )
    if (allocated(this%qm)) deallocate(this%qm)

    this%is_initialized = .false.

  end subroutine tracers_clear

  subroutine tracers_final(this)

    type(tracers_type), intent(inout) :: this

    call this%clear()

  end subroutine tracers_final

end module tracer_types_mod