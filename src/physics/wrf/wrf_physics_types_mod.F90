module wrf_physics_types_mod

  use const_mod
  use physics_types_mod
  use wrf_namelist_mod
  use lsm_noahmp_types_mod

  implicit none

  private

  public wrf_state_type

  type, extends(physics_state_type) :: wrf_state_type
    ! Time step index
    integer :: time_step = 0
    ! U wind component at 10 m
    real(r8), allocatable, dimension(:) :: u10
    ! V wind component at 10 m
    real(r8), allocatable, dimension(:) :: v10
    ! Potential temperature at 2 m
    real(r8), allocatable, dimension(:) :: pt2
    ! Specific humidity at 2 m
    real(r8), allocatable, dimension(:) :: qv2
    type(lsm_noahmp_state_type) lsm_noahmp
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

    allocate(this%u10(mesh%ncol))
    allocate(this%v10(mesh%ncol))
    allocate(this%pt2(mesh%ncol))
    allocate(this%qv2(mesh%ncol))

    call this%lsm_noahmp%init(mesh)

    call this%physics_state_init(mesh)

  end subroutine wrf_state_init

  subroutine wrf_state_clear(this)

    class(wrf_state_type), intent(inout) :: this

    this%time_step = 0

    if (allocated(this%u10)) deallocate(this%u10)
    if (allocated(this%v10)) deallocate(this%v10)
    if (allocated(this%pt2)) deallocate(this%pt2)
    if (allocated(this%qv2)) deallocate(this%qv2)

    call this%lsm_noahmp%clear()

    call this%physics_state_clear()

  end subroutine wrf_state_clear

  subroutine wrf_state_final(this)

    type(wrf_state_type), intent(inout) :: this

    call this%clear()

  end subroutine wrf_state_final

end module wrf_physics_types_mod
