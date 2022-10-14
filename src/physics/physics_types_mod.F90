module physics_types_mod

  use const_mod
  use namelist_mod
  use mesh_mod

  implicit none

  private

  public pstate_type
  public ptend_type

  type pstate_type
    integer :: ncol   = 0
    integer :: nlev   = 0
    integer , allocatable, dimension(:  ) :: i
    integer , allocatable, dimension(:  ) :: j
    real(r8), allocatable, dimension(:  ) :: lon
    real(r8), allocatable, dimension(:  ) :: lat
    real(r8), allocatable, dimension(:,:) :: u
    real(r8), allocatable, dimension(:,:) :: v
    real(r8), allocatable, dimension(:,:) :: t
    real(r8), allocatable, dimension(:,:) :: pt     ! Potential temperature
    real(r8), allocatable, dimension(:,:) :: p
    real(r8), allocatable, dimension(:,:) :: dp
    real(r8), allocatable, dimension(:,:) :: rdp    ! 1 / dp
    real(r8), allocatable, dimension(:,:) :: p_lev
    real(r8), allocatable, dimension(:,:) :: z
    real(r8), allocatable, dimension(:,:) :: z_lev
    real(r8), allocatable, dimension(:,:) :: sh     ! Specific humidity
    real(r8), allocatable, dimension(:,:) :: rho    ! Air density
    real(r8), allocatable, dimension(:  ) :: ps
    real(r8), allocatable, dimension(:  ) :: precl  ! Large scale precipitation
  contains
    procedure :: init  => pstate_init
    procedure :: clear => pstate_clear
    final :: pstate_final
  end type pstate_type

  type ptend_type
    integer :: ncol = 0
    integer :: nlev = 0
    real(r8), allocatable, dimension(:,:) :: dudt
    real(r8), allocatable, dimension(:,:) :: dvdt
    real(r8), allocatable, dimension(:,:) :: dtdt
    real(r8), allocatable, dimension(:,:) :: dshdt
    logical :: updated_u  = .false.
    logical :: updated_v  = .false.
    logical :: updated_t  = .false.
    logical :: updated_sh = .false.
  contains
    procedure :: init  => ptend_init
    procedure :: clear => ptend_clear
    procedure :: reset => ptend_reset
    final :: ptend_final
  end type ptend_type

contains

  subroutine pstate_init(this, mesh)

    class(pstate_type), intent(inout) :: this
    type(mesh_type), intent(in) :: mesh

    integer i, j, icol

    call this%clear()

    this%ncol = mesh%num_full_lon * mesh%num_full_lat
    this%nlev = mesh%num_full_lev
    allocate(this%i     (this%ncol            ))
    allocate(this%j     (this%ncol            ))
    allocate(this%lon   (this%ncol            ))
    allocate(this%lat   (this%ncol            ))
    allocate(this%p_lev (this%ncol,this%nlev+1))
    allocate(this%z_lev (this%ncol,this%nlev+1))
    allocate(this%u     (this%ncol,this%nlev  ))
    allocate(this%v     (this%ncol,this%nlev  ))
    allocate(this%t     (this%ncol,this%nlev  ))
    allocate(this%pt    (this%ncol,this%nlev  ))
    allocate(this%p     (this%ncol,this%nlev  ))
    allocate(this%dp    (this%ncol,this%nlev  ))
    allocate(this%rdp   (this%ncol,this%nlev  ))
    allocate(this%z     (this%ncol,this%nlev  ))
    allocate(this%sh    (this%ncol,this%nlev  ))
    allocate(this%rho   (this%ncol,this%nlev  ))
    allocate(this%ps    (this%ncol            ))
    allocate(this%precl (this%ncol            ))

    icol = 0
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        icol = icol + 1
        this%i  (icol) = i
        this%j  (icol) = j
        this%lon(icol) = mesh%full_lon(i)
        this%lat(icol) = mesh%full_lat(j)
      end do
    end do

  end subroutine pstate_init

  subroutine pstate_clear(this)

    class(pstate_type), intent(inout) :: this

    if (allocated(this%i     )) deallocate(this%i    )
    if (allocated(this%j     )) deallocate(this%j    )
    if (allocated(this%lon   )) deallocate(this%lon  )
    if (allocated(this%lat   )) deallocate(this%lat  )
    if (allocated(this%u     )) deallocate(this%u    )
    if (allocated(this%v     )) deallocate(this%v    )
    if (allocated(this%t     )) deallocate(this%t    )
    if (allocated(this%pt    )) deallocate(this%pt   )
    if (allocated(this%p     )) deallocate(this%p    )
    if (allocated(this%dp    )) deallocate(this%dp   )
    if (allocated(this%rdp   )) deallocate(this%rdp  )
    if (allocated(this%p_lev )) deallocate(this%p_lev)
    if (allocated(this%z     )) deallocate(this%z    )
    if (allocated(this%z_lev )) deallocate(this%z_lev)
    if (allocated(this%sh    )) deallocate(this%sh   )
    if (allocated(this%rho   )) deallocate(this%rho  )
    if (allocated(this%ps    )) deallocate(this%ps   )
    if (allocated(this%precl )) deallocate(this%precl)

  end subroutine pstate_clear

  subroutine pstate_final(this)

    type(pstate_type), intent(inout) :: this

    call this%clear()

  end subroutine pstate_final

  subroutine ptend_init(this, mesh)

    class(ptend_type), intent(inout) :: this
    type(mesh_type), intent(in) :: mesh

    call this%clear()

    this%ncol = mesh%num_full_lon * mesh%num_full_lat
    this%nlev = mesh%num_full_lev
    allocate(this%dvdt (this%ncol,this%nlev))
    allocate(this%dtdt (this%ncol,this%nlev))
    allocate(this%dshdt(this%ncol,this%nlev))
    allocate(this%dudt (this%ncol,this%nlev))    

  end subroutine ptend_init

  subroutine ptend_clear(this)

    class(ptend_type), intent(inout) :: this

    if (allocated(this%dudt )) deallocate(this%dudt )
    if (allocated(this%dvdt )) deallocate(this%dvdt )
    if (allocated(this%dtdt )) deallocate(this%dtdt )
    if (allocated(this%dshdt)) deallocate(this%dshdt)

  end subroutine ptend_clear

  subroutine ptend_final(this)

    type(ptend_type), intent(inout) :: this

    call this%clear()

  end subroutine ptend_final

  subroutine ptend_reset(this)

    class(ptend_type), intent(inout) :: this

    this%dudt  = 0; this%updated_u  = .false.
    this%dvdt  = 0; this%updated_v  = .false.
    this%dtdt  = 0; this%updated_t  = .false.
    this%dshdt = 0; this%updated_sh = .false.

  end subroutine ptend_reset

end module physics_types_mod
