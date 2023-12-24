! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module physics_types_mod

  use const_mod
  use namelist_mod
  use physics_mesh_mod
  use tracer_types_mod

  implicit none

  private

  public physics_mesh_type
  public physics_state_type
  public physics_tend_type
  public physics_static_type

  type, abstract :: physics_state_type
    type(physics_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:  ) :: u      ! U-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: v      ! V-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: t      ! Temperature (K)
    real(r8), allocatable, dimension(:,:,:) :: q      ! Tracer mixing ratio (moist)
    real(r8), allocatable, dimension(:,:  ) :: p      ! Full pressure (hydrostatic) on full level (Pa)
    real(r8), allocatable, dimension(:,:  ) :: p_lev  ! Full pressure (hydrostatic) on half level (Pa)
    real(r8), allocatable, dimension(:,:  ) :: pk     ! Exner function of full pressure (hydrostatic) on full level
    real(r8), allocatable, dimension(:,:  ) :: pk_lev ! Exner function of full pressure (hydrostatic) on half levels
    real(r8), allocatable, dimension(:,:  ) :: dp     ! Full pressure thickness (Pa)
    real(r8), allocatable, dimension(:,:  ) :: omg    ! Vertical pressure velocity (Pa s-1)
    real(r8), allocatable, dimension(:,:  ) :: z      ! Height on full levels
    real(r8), allocatable, dimension(:,:  ) :: z_lev  ! Height on half levels
    real(r8), allocatable, dimension(:,:  ) :: dz     ! Height thickness on full levels
    real(r8), allocatable, dimension(:,:  ) :: rho    ! Air density
    real(r8), allocatable, dimension(:    ) :: ps     ! Surface pressure (Pa)
    real(r8), allocatable, dimension(:    ) :: ts     ! Surface temperature (K)
    real(r8), allocatable, dimension(:    ) :: wsb    ! Wind speed on lowest model level (m s-1)
    real(r8), allocatable, dimension(:    ) :: ustar  ! u* in similarity theory (m s-1)
  contains
    procedure physics_state_init
    procedure physics_state_clear
  end type physics_state_type

  type, abstract :: physics_tend_type
    type(physics_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:  ) :: dudt
    real(r8), allocatable, dimension(:,:  ) :: dvdt
    real(r8), allocatable, dimension(:,:  ) :: dtdt
    real(r8), allocatable, dimension(:,:  ) :: dptdt
    real(r8), allocatable, dimension(:,:,:) :: dqdt
    logical :: updated_u  = .false.
    logical :: updated_v  = .false.
    logical :: updated_t  = .false.
    logical, allocatable :: updated_q(:)
  contains
    procedure physics_tend_init
    procedure physics_tend_clear
    procedure physics_tend_reset
  end type physics_tend_type

  type, abstract :: physics_static_type
    type(physics_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:) :: z0      ! Surface roughness length (m)
    real(r8), allocatable, dimension(:) :: sfc_alb ! Surface albedo
  contains
    procedure physics_static_init
    procedure physics_static_clear
  end type physics_static_type

contains

  subroutine physics_state_init(this, mesh)

    class(physics_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%physics_state_clear()

    this%mesh => mesh

    allocate(this%u         (mesh%ncol,mesh%nlev         ))
    allocate(this%v         (mesh%ncol,mesh%nlev         ))
    allocate(this%t         (mesh%ncol,mesh%nlev         ))
    allocate(this%q         (mesh%ncol,mesh%nlev,ntracers))
    allocate(this%p         (mesh%ncol,mesh%nlev         ))
    allocate(this%p_lev     (mesh%ncol,mesh%nlev+1       ))
    allocate(this%pk        (mesh%ncol,mesh%nlev         ))
    allocate(this%pk_lev    (mesh%ncol,mesh%nlev+1       ))
    allocate(this%dp        (mesh%ncol,mesh%nlev         ))
    allocate(this%omg       (mesh%ncol,mesh%nlev         ))
    allocate(this%z         (mesh%ncol,mesh%nlev         ))
    allocate(this%z_lev     (mesh%ncol,mesh%nlev+1       ))
    allocate(this%dz        (mesh%ncol,mesh%nlev         ))
    allocate(this%rho       (mesh%ncol,mesh%nlev         ))
    allocate(this%ps        (mesh%ncol                   ))
    allocate(this%ts        (mesh%ncol                   ))
    allocate(this%wsb       (mesh%ncol                   ))
    allocate(this%ustar     (mesh%ncol                   ))

  end subroutine physics_state_init

  subroutine physics_state_clear(this)

    class(physics_state_type), intent(inout) :: this

    this%mesh => null()

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%t        )) deallocate(this%t        )
    if (allocated(this%q        )) deallocate(this%q        )
    if (allocated(this%p        )) deallocate(this%p        )
    if (allocated(this%p_lev    )) deallocate(this%p_lev    )
    if (allocated(this%pk       )) deallocate(this%pk       )
    if (allocated(this%pk_lev   )) deallocate(this%pk_lev   )
    if (allocated(this%dp       )) deallocate(this%dp       )
    if (allocated(this%omg      )) deallocate(this%omg      )
    if (allocated(this%z        )) deallocate(this%z        )
    if (allocated(this%z_lev    )) deallocate(this%z_lev    )
    if (allocated(this%dz       )) deallocate(this%dz       )
    if (allocated(this%rho      )) deallocate(this%rho      )
    if (allocated(this%ps       )) deallocate(this%ps       )
    if (allocated(this%ts       )) deallocate(this%ts       )
    if (allocated(this%wsb      )) deallocate(this%wsb      )
    if (allocated(this%ustar    )) deallocate(this%ustar    )

  end subroutine physics_state_clear

  subroutine physics_tend_init(this, mesh)

    class(physics_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%physics_tend_clear()

    this%mesh => mesh

    allocate(this%dudt (mesh%ncol,mesh%nlev         ))
    allocate(this%dvdt (mesh%ncol,mesh%nlev         ))
    allocate(this%dtdt (mesh%ncol,mesh%nlev         ))
    allocate(this%dptdt(mesh%ncol,mesh%nlev         ))
    allocate(this%dqdt (mesh%ncol,mesh%nlev,ntracers))
    allocate(this%updated_q(ntracers))

  end subroutine physics_tend_init

  subroutine physics_tend_clear(this)

    class(physics_tend_type), intent(inout) :: this

    this%mesh => null()

    if (allocated(this%dudt )) deallocate(this%dudt )
    if (allocated(this%dvdt )) deallocate(this%dvdt )
    if (allocated(this%dtdt )) deallocate(this%dtdt )
    if (allocated(this%dptdt)) deallocate(this%dptdt)
    if (allocated(this%dqdt )) deallocate(this%dqdt )
    if (allocated(this%updated_q)) deallocate(this%updated_q)

  end subroutine physics_tend_clear

  subroutine physics_tend_reset(this)

    class(physics_tend_type), intent(inout) :: this

    this%dudt  = 0; this%updated_u  = .false.
    this%dvdt  = 0; this%updated_v  = .false.
    this%dtdt  = 0; this%updated_t  = .false.
    this%dqdt  = 0; this%updated_q  = .false.

  end subroutine physics_tend_reset

  subroutine physics_static_init(this, mesh)

    class(physics_static_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%physics_static_clear()

    this%mesh => mesh

    allocate(this%z0     (mesh%ncol))
    allocate(this%sfc_alb(mesh%ncol))

  end subroutine physics_static_init

  subroutine physics_static_clear(this)

    class(physics_static_type), intent(inout) :: this

    this%mesh => null()
    if (allocated(this%z0     )) deallocate(this%z0     )
    if (allocated(this%sfc_alb)) deallocate(this%sfc_alb)

  end subroutine physics_static_clear

end module physics_types_mod
