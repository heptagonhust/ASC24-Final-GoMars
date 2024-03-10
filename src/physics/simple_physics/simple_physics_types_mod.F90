! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module simple_physics_types_mod

  use const_mod
  use tracer_mod
  use physics_types_mod

  implicit none

  type, extends(physics_state_type) :: simple_state_type
    real(r8), pointer    , dimension(:,:) :: qv
    real(r8), allocatable, dimension(:  ) :: precl
  contains
    procedure :: init  => simple_state_init
    procedure :: clear => simple_state_clear
    final simple_state_final
  end type simple_state_type

  type, extends(physics_tend_type) :: simple_tend_type
    real(r8), pointer    , dimension(:,:) :: dqvdt
  contains
    procedure :: init  => simple_tend_init
    procedure :: clear => simple_tend_clear
    final simple_tend_final
  end type simple_tend_type

contains

  subroutine simple_state_init(this, mesh)

    class(simple_state_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_state_init(mesh)

    this%qv => this%q(:,:,idx_qv)
    allocate(this%precl(mesh%ncol))

  end subroutine simple_state_init

  subroutine simple_state_clear(this)

    class(simple_state_type), intent(inout) :: this

    call this%physics_state_clear()

    if (allocated(this%precl)) deallocate(this%precl)

  end subroutine simple_state_clear

  subroutine simple_state_final(this)

    type(simple_state_type), intent(inout) :: this

    call this%clear()

  end subroutine simple_state_final

  subroutine simple_tend_init(this, mesh)

    class(simple_tend_type), intent(inout), target :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

    this%dqvdt => this%dqdt(:,:,idx_qv)

  end subroutine simple_tend_init

  subroutine simple_tend_clear(this)

    class(simple_tend_type), intent(inout) :: this

    call this%physics_tend_clear()

  end subroutine simple_tend_clear

  subroutine simple_tend_final(this)

    type(simple_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine simple_tend_final

end module simple_physics_types_mod
