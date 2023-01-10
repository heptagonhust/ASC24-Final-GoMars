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
    real(r8), allocatable, dimension(:,:) :: p_exn  ! Exner function
    real(r8), allocatable, dimension(:,:) :: dp     ! Pressure difference
    real(r8), allocatable, dimension(:,:) :: rdp    ! 1 / dp
    real(r8), allocatable, dimension(:,:) :: p_lev
    real(r8), allocatable, dimension(:,:) :: z
    real(r8), allocatable, dimension(:,:) :: dz
    real(r8), allocatable, dimension(:,:) :: z_lev
    real(r8), allocatable, dimension(:,:) :: sh     ! Specific humidity
    real(r8), allocatable, dimension(:,:) :: qv     ! Water vapor mixing ratio
    real(r8), allocatable, dimension(:,:) :: qc     ! Cloud water mixing ratio
    real(r8), allocatable, dimension(:,:) :: qi     ! Cloud ice mixing ratio
    real(r8), allocatable, dimension(:,:) :: rho    ! Air density
    real(r8), allocatable, dimension(:  ) :: emis   ! Surface emissivity
    real(r8), allocatable, dimension(:  ) :: alb    ! Surface albedo
    real(r8), allocatable, dimension(:  ) :: ps     ! Surface pressure (Pa)
    real(r8), allocatable, dimension(:  ) :: ts     ! Surface temperature (K)
    real(r8), allocatable, dimension(:  ) :: precl  ! Large scale precipitation
    real(r8), allocatable, dimension(:  ) :: pblh   ! PBL height (m)
    integer , allocatable, dimension(:  ) :: pblk   ! PBL level index
    real(r8), allocatable, dimension(:  ) :: N      ! Brunt-Väisälä frequency (s-1)
    real(r8), allocatable, dimension(:  ) :: z0     ! Roughness height
    real(r8), allocatable, dimension(:  ) :: ustar  ! u* in similarity theory (m s-1)
    real(r8), allocatable, dimension(:  ) :: wstar  ! Mixed-layer velocity scale (m s-1)
    real(r8), allocatable, dimension(:  ) :: wsp    ! Wind speed at lowest model level (m s-1)
    real(r8), allocatable, dimension(:  ) :: u10    ! U-wind speed at 10m (m s-1)
    real(r8), allocatable, dimension(:  ) :: v10    ! V-wind speed at 10m (m s-1)
    real(r8), allocatable, dimension(:  ) :: uos    ! Sea surface zonal current (m s-1)
    real(r8), allocatable, dimension(:  ) :: vos    ! Sea surface meridional current (m s-1)
    real(r8), allocatable, dimension(:  ) :: Rib    ! Bulk Richardson number in surface layer
    real(r8), allocatable, dimension(:  ) :: psim   ! Similarity stability function for momentum
    real(r8), allocatable, dimension(:  ) :: psih   ! Similarity stability function for heat
    real(r8), allocatable, dimension(:  ) :: land   ! Land mask (1 for land, 2 for water)
    real(r8), allocatable, dimension(:  ) :: hfx    ! Upward heat flux at surface (W m-2)
    real(r8), allocatable, dimension(:  ) :: qfx    ! Upward moisture flux at surface (kg s-1 m-2)
    real(r8), allocatable, dimension(:,:) :: exch_h ! Exchange coefficient for heat (K m s-1)
    real(r8), allocatable, dimension(:,:) :: delta  ! Entrainment layer depth (m)
    real(r8), allocatable, dimension(:  ) :: co2ice ! CO2 ice on the surface (Kg m-2)
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
    real(r8), allocatable, dimension(:,:) :: dptdt
    real(r8), allocatable, dimension(:,:) :: dptdt_rad
    real(r8), allocatable, dimension(:,:) :: dshdt
    real(r8), allocatable, dimension(:,:) :: dqvdt
    real(r8), allocatable, dimension(:,:) :: dqcdt
    real(r8), allocatable, dimension(:,:) :: dqidt
    logical :: updated_u  = .false.
    logical :: updated_v  = .false.
    logical :: updated_t  = .false.
    logical :: updated_pt = .false.
    logical :: updated_sh = .false.
    logical :: updated_qv = .false.
    logical :: updated_qc = .false.
    logical :: updated_qi = .false.
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

    this%ncol = mesh%full_nlon * mesh%full_nlat
    this%nlev = mesh%full_nlev
    allocate(this%i         (this%ncol            ))
    allocate(this%j         (this%ncol            ))
    allocate(this%lon       (this%ncol            ))
    allocate(this%lat       (this%ncol            ))
    allocate(this%u         (this%ncol,this%nlev  ))
    allocate(this%v         (this%ncol,this%nlev  ))
    allocate(this%t         (this%ncol,this%nlev  ))
    allocate(this%pt        (this%ncol,this%nlev  ))
    allocate(this%p         (this%ncol,this%nlev  ))
    allocate(this%p_exn     (this%ncol,this%nlev  ))
    allocate(this%dp        (this%ncol,this%nlev  ))
    allocate(this%rdp       (this%ncol,this%nlev  ))
    allocate(this%p_lev     (this%ncol,this%nlev+1))
    allocate(this%z         (this%ncol,this%nlev  ))
    allocate(this%dz        (this%ncol,this%nlev  ))
    allocate(this%z_lev     (this%ncol,this%nlev+1))
    allocate(this%sh        (this%ncol,this%nlev  ))
    allocate(this%qv        (this%ncol,this%nlev  ))
    allocate(this%qc        (this%ncol,this%nlev  ))
    allocate(this%qi        (this%ncol,this%nlev  ))
    allocate(this%rho       (this%ncol,this%nlev  ))
    allocate(this%emis      (this%ncol            ))
    allocate(this%alb       (this%ncol            ))
    allocate(this%ps        (this%ncol            ))
    allocate(this%ts        (this%ncol            ))
    allocate(this%precl     (this%ncol            ))
    allocate(this%pblh      (this%ncol            ))
    allocate(this%pblk      (this%ncol            ))
    allocate(this%N         (this%ncol            ))
    allocate(this%z0        (this%ncol            ))
    allocate(this%ustar     (this%ncol            ))
    allocate(this%wstar     (this%ncol            ))
    allocate(this%wsp       (this%ncol            ))
    allocate(this%u10       (this%ncol            ))
    allocate(this%v10       (this%ncol            ))
    allocate(this%uos       (this%ncol            ))
    allocate(this%vos       (this%ncol            ))
    allocate(this%Rib       (this%ncol            ))
    allocate(this%psim      (this%ncol            ))
    allocate(this%psih      (this%ncol            ))
    allocate(this%land      (this%ncol            ))
    allocate(this%hfx       (this%ncol            ))
    allocate(this%qfx       (this%ncol            ))
    allocate(this%exch_h    (this%ncol,this%nlev  ))
    allocate(this%delta     (this%ncol,this%nlev  ))

    select case (planet)
    case ('mars')
      allocate(this%co2ice  (this%ncol            ))
    end select

    icol = 0
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
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

    if (allocated(this%i        )) deallocate(this%i        )
    if (allocated(this%j        )) deallocate(this%j        )
    if (allocated(this%lon      )) deallocate(this%lon      )
    if (allocated(this%lat      )) deallocate(this%lat      )
    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%t        )) deallocate(this%t        )
    if (allocated(this%pt       )) deallocate(this%pt       )
    if (allocated(this%p        )) deallocate(this%p        )
    if (allocated(this%p_exn    )) deallocate(this%p_exn    )
    if (allocated(this%dp       )) deallocate(this%dp       )
    if (allocated(this%rdp      )) deallocate(this%rdp      )
    if (allocated(this%p_lev    )) deallocate(this%p_lev    )
    if (allocated(this%z        )) deallocate(this%z        )
    if (allocated(this%dz       )) deallocate(this%dz       )
    if (allocated(this%z_lev    )) deallocate(this%z_lev    )
    if (allocated(this%sh       )) deallocate(this%sh       )
    if (allocated(this%qv       )) deallocate(this%qv       )
    if (allocated(this%qc       )) deallocate(this%qc       )
    if (allocated(this%qi       )) deallocate(this%qi       )
    if (allocated(this%rho      )) deallocate(this%rho      )
    if (allocated(this%ps       )) deallocate(this%ps       )
    if (allocated(this%precl    )) deallocate(this%precl    )
    if (allocated(this%pblh     )) deallocate(this%pblh     )
    if (allocated(this%pblk     )) deallocate(this%pblk     )
    if (allocated(this%N        )) deallocate(this%N        )
    if (allocated(this%z0       )) deallocate(this%z0       )
    if (allocated(this%ustar    )) deallocate(this%ustar    )
    if (allocated(this%wstar    )) deallocate(this%wstar    )
    if (allocated(this%wsp      )) deallocate(this%wsp      )
    if (allocated(this%u10      )) deallocate(this%u10      )
    if (allocated(this%v10      )) deallocate(this%v10      )
    if (allocated(this%uos      )) deallocate(this%uos      )
    if (allocated(this%vos      )) deallocate(this%vos      )
    if (allocated(this%Rib      )) deallocate(this%Rib      )
    if (allocated(this%psim     )) deallocate(this%psim     )
    if (allocated(this%psih     )) deallocate(this%psih     )
    if (allocated(this%land     )) deallocate(this%land     )
    if (allocated(this%hfx      )) deallocate(this%hfx      )
    if (allocated(this%qfx      )) deallocate(this%qfx      )
    if (allocated(this%exch_h   )) deallocate(this%exch_h   )
    if (allocated(this%delta    )) deallocate(this%delta    )

  end subroutine pstate_clear

  subroutine pstate_final(this)

    type(pstate_type), intent(inout) :: this

    call this%clear()

  end subroutine pstate_final

  subroutine ptend_init(this, mesh)

    class(ptend_type), intent(inout) :: this
    type(mesh_type), intent(in) :: mesh

    call this%clear()

    this%ncol = mesh%full_nlon * mesh%full_nlat
    this%nlev = mesh%full_nlev
    allocate(this%dudt      (this%ncol,this%nlev))
    allocate(this%dvdt      (this%ncol,this%nlev))
    allocate(this%dtdt      (this%ncol,this%nlev))
    allocate(this%dptdt     (this%ncol,this%nlev))
    allocate(this%dptdt_rad (this%ncol,this%nlev)); this%dptdt_rad = 0
    allocate(this%dshdt     (this%ncol,this%nlev))
    allocate(this%dqvdt     (this%ncol,this%nlev))
    allocate(this%dqcdt     (this%ncol,this%nlev))
    allocate(this%dqidt     (this%ncol,this%nlev))

  end subroutine ptend_init

  subroutine ptend_clear(this)

    class(ptend_type), intent(inout) :: this

    if (allocated(this%dudt     )) deallocate(this%dudt     )
    if (allocated(this%dvdt     )) deallocate(this%dvdt     )
    if (allocated(this%dtdt     )) deallocate(this%dtdt     )
    if (allocated(this%dptdt    )) deallocate(this%dptdt    )
    if (allocated(this%dptdt_rad)) deallocate(this%dptdt_rad)
    if (allocated(this%dshdt    )) deallocate(this%dshdt    )
    if (allocated(this%dqvdt    )) deallocate(this%dqvdt    )
    if (allocated(this%dqcdt    )) deallocate(this%dqcdt    )
    if (allocated(this%dqidt    )) deallocate(this%dqidt    )

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
    this%dptdt = 0; this%updated_pt = .false.
    this%dshdt = 0; this%updated_sh = .false.
    this%dqvdt = 0; this%updated_qv = .false.
    this%dqcdt = 0; this%updated_qc = .false.
    this%dqidt = 0; this%updated_qi = .false.

  end subroutine ptend_reset

end module physics_types_mod
