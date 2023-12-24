module physics_mesh_mod

  use const_mod

  implicit none

  private

  public physics_mesh_type

  type physics_mesh_type
    integer :: ncol = 0
    integer :: nlev = 0
    real(r8), allocatable, dimension(:) :: lon
    real(r8), allocatable, dimension(:) :: lat
    real(r8), allocatable, dimension(:) :: area
  contains
    procedure :: init => physics_mesh_init
    procedure :: clear => physics_mesh_clear
    final physics_mesh_final
  end type physics_mesh_type

contains

  subroutine physics_mesh_init(this, ncol, nlev, lon, lat, area)

    class(physics_mesh_type), intent(inout) :: this
    integer , intent(in) :: ncol
    integer , intent(in) :: nlev
    real(r8), intent(in) :: lon (ncol)
    real(r8), intent(in) :: lat (ncol)
    real(r8), intent(in) :: area(ncol)

    call this%clear()

    this%ncol = ncol
    this%nlev = nlev
    allocate(this%lon (ncol)); this%lon  = lon
    allocate(this%lat (ncol)); this%lat  = lat
    allocate(this%area(ncol)); this%area = area

  end subroutine physics_mesh_init

  subroutine physics_mesh_clear(this)

    class(physics_mesh_type), intent(inout) :: this

    this%ncol = 0
    this%nlev = 0
    if (allocated(this%lon )) deallocate(this%lon )
    if (allocated(this%lat )) deallocate(this%lat )
    if (allocated(this%area)) deallocate(this%area)

  end subroutine physics_mesh_clear

  subroutine physics_mesh_final(this)

    type(physics_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine physics_mesh_final

end module physics_mesh_mod
