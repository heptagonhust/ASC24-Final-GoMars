module physics_types_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use tracer_types_mod

  implicit none

  private

  public pstate_type
  public ptend_type

  type pstate_type
    integer :: ncol = 0
    integer :: nlev = 0
    integer , allocatable, dimension(:    ) :: i
    integer , allocatable, dimension(:    ) :: j
    real(r8), allocatable, dimension(:    ) :: lon
    real(r8), allocatable, dimension(:    ) :: lat
    real(r8), allocatable, dimension(:    ) :: area     ! Cell area (m2)
    ! Wind
    real(r8), allocatable, dimension(:,:  ) :: u        ! U-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: u_new    ! Updated u-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: v        ! V-wind speed (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: v_new    ! Updated v-wind speed (m s-1)
    real(r8), allocatable, dimension(:    ) :: wsb      ! Wind speed on lowest model level (m s-1)
    real(r8), allocatable, dimension(:    ) :: u10      ! U-wind speed at 10m (m s-1)
    real(r8), allocatable, dimension(:    ) :: v10      ! V-wind speed at 10m (m s-1)
    ! Temperature
    real(r8), allocatable, dimension(:,:  ) :: t        ! Temperature (K)
    real(r8), allocatable, dimension(:,:  ) :: t_new    ! Temperature (K)
    real(r8), allocatable, dimension(:,:  ) :: tv       ! Virtual temperature (K)
    real(r8), allocatable, dimension(:,:  ) :: pt       ! Potential temperature (K)
    real(r8), allocatable, dimension(:,:  ) :: ptv      ! Virtual potential temperature (K)
    real(r8), allocatable, dimension(:    ) :: t_sfc    ! Surface or ground temperature (K)
    ! Pressure
    real(r8), allocatable, dimension(:,:  ) :: p        ! Full pressure (hydrostatic) on full levels (Pa)
    real(r8), allocatable, dimension(:,:  ) :: p_lev    ! Full pressure (hydrostatic) on half levels (Pa)
    real(r8), allocatable, dimension(:,:  ) :: pk       ! Exner function of full pressure (hydrostatic) on full levels
    real(r8), allocatable, dimension(:,:  ) :: pk_lev   ! Exner function of full pressure (hydrostatic) on half levels
    real(r8), allocatable, dimension(:,:  ) :: dp       ! Full pressure thickness (Pa)
    real(r8), allocatable, dimension(:,:  ) :: rdp      ! 1 / dp
    real(r8), allocatable, dimension(:,:  ) :: lnp_lev  ! Logrithm of full pressure on half levels
    real(r8), allocatable, dimension(:,:  ) :: omg      ! Vertical pressure velocity (Pa s-1)
    ! Height
    real(r8), allocatable, dimension(:,:  ) :: z        ! Height on full levels
    real(r8), allocatable, dimension(:,:  ) :: z_lev    ! Height on half levels
    real(r8), allocatable, dimension(:,:  ) :: dz       ! Height thickness on full levels
    ! Tracers
    ! NOTE: Dynamical core uses dry mixing ratios, but physics uses moist mixing ratios.
    real(r8), allocatable, dimension(:,:,:) :: q        ! Tracer mixing ratio (moist)
    ! Moisture
    real(r8), pointer    , dimension(:,:  ) :: qv       ! Water vapor mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qc       ! Cloud water mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qi       ! Cloud ice mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qr       ! Rain mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qs       ! Snow mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qg       ! Grauple mixing ratio (moist)
    real(r8), pointer    , dimension(:,:  ) :: qh       ! Hail mixing ratio (moist)
    ! Ozone
    real(r8), pointer    , dimension(:,:  ) :: qo3      ! Ozone mixing ratio (moist)
    ! Stability
    real(r8), allocatable, dimension(:,:  ) :: n2_lev   ! Square of Brunt-Väisälä frequency (s-2) on half levels
    real(r8), allocatable, dimension(:,:  ) :: ri_lev   ! Local Richardson number on half levels
    ! Surface layer
    real(r8), allocatable, dimension(:    ) :: emis     ! Surface emissivity
    real(r8), allocatable, dimension(:    ) :: alb      ! Surface albedo
    real(r8), allocatable, dimension(:    ) :: ps       ! Surface pressure (Pa)
    real(r8), allocatable, dimension(:    ) :: ts       ! Surface temperature (K)
    real(r8), allocatable, dimension(:    ) :: land     ! Land mask (1 for land, 2 for water)
    real(r8), allocatable, dimension(:    ) :: hfx      ! Upward heat flux at surface (W m-2)
    real(r8), allocatable, dimension(:    ) :: qfx      ! Upward moisture flux at surface (kg s-1 m-2)
    real(r8), allocatable, dimension(:    ) :: z0       ! Roughness height
    real(r8), allocatable, dimension(:    ) :: ustar    ! u* in similarity theory (m s-1)
    real(r8), allocatable, dimension(:    ) :: ptstar   ! pt* (K)
    real(r8), allocatable, dimension(:    ) :: psim     ! Similarity stability function for momentum
    real(r8), allocatable, dimension(:    ) :: psih     ! Similarity stability function for heat
    real(r8), allocatable, dimension(:    ) :: rib      ! Bulk Richardson number in surface layer
    real(r8), allocatable, dimension(:    ) :: uos      ! Sea surface zonal current (m s-1)
    real(r8), allocatable, dimension(:    ) :: vos      ! Sea surface meridional current (m s-1)
    ! Boundary layer
    real(r8), allocatable, dimension(:    ) :: pblh     ! PBL height (m)
    integer , allocatable, dimension(:    ) :: pblk     ! PBL level index
    real(r8), allocatable, dimension(:    ) :: wstar    ! Mixed-layer velocity scale (m s-1)
    real(r8), allocatable, dimension(:,:  ) :: exch_h   ! Exchange coefficient for heat (K m s-1)
    real(r8), allocatable, dimension(:,:  ) :: delta    ! Entrainment layer depth (m)
    ! Precipitation
    real(r8), allocatable, dimension(:    ) :: precl    ! Large scale precipitation
    ! Mars
    real(r8), allocatable, dimension(:    ) :: co2ice   ! CO2 ice on the surface (Kg m-2)
    !
    real(r8), allocatable, dimension(:,:  ) :: rho      ! Air density
    real(r8), allocatable, dimension(:,:  ) :: cp       ! Specific heat capacity of total air in constant pressure
    real(r8), allocatable, dimension(:,:  ) :: cv       ! Specific heat capacity of total air in constant volume
    real(r8), allocatable, dimension(:,:  ) :: tep      ! Total enery in cp * T + gz + K
    real(r8), allocatable, dimension(:,:  ) :: tev      ! Total enery in cv * T + gz + K
  contains
    procedure :: init => pstate_init
    procedure :: clear => pstate_clear
    final :: pstate_final
  end type pstate_type

  type ptend_type
    integer :: ncol = 0
    integer :: nlev = 0
    real(r8), allocatable, dimension(:,:  ) :: dudt
    real(r8), allocatable, dimension(:,:  ) :: dvdt
    real(r8), allocatable, dimension(:,:  ) :: dtdt
    real(r8), allocatable, dimension(:,:,:) :: dqdt
    real(r8), allocatable, dimension(:,:  ) :: dptdt
    real(r8), allocatable, dimension(:,:  ) :: dptdt_rad
    logical :: updated_u  = .false.
    logical :: updated_v  = .false.
    logical :: updated_t  = .false.
    logical :: updated_pt = .false.
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

    class(pstate_type), intent(inout), target :: this
    type(mesh_type), intent(in) :: mesh

    integer i, j, icol

    call this%clear()

    this%ncol = mesh%full_nlon * mesh%full_nlat
    this%nlev = mesh%full_nlev
    allocate(this%i         (this%ncol            ))
    allocate(this%j         (this%ncol            ))
    allocate(this%lon       (this%ncol            ))
    allocate(this%lat       (this%ncol            ))
    allocate(this%area      (this%ncol            ))
    ! Wind
    allocate(this%u         (this%ncol,this%nlev  ))
    allocate(this%v         (this%ncol,this%nlev  ))
    allocate(this%wsb       (this%ncol            ))
    allocate(this%u10       (this%ncol            ))
    allocate(this%v10       (this%ncol            ))
    ! Temperature
    allocate(this%t         (this%ncol,this%nlev  ))
    allocate(this%tv        (this%ncol,this%nlev  ))
    allocate(this%pt        (this%ncol,this%nlev  ))
    allocate(this%ptv       (this%ncol,this%nlev  ))
    allocate(this%t_sfc     (this%ncol            ))
    ! Pressure
    allocate(this%p         (this%ncol,this%nlev  ))
    allocate(this%p_lev     (this%ncol,this%nlev+1))
    allocate(this%pk        (this%ncol,this%nlev  ))
    allocate(this%pk_lev    (this%ncol,this%nlev+1))
    allocate(this%dp        (this%ncol,this%nlev  ))
    allocate(this%rdp       (this%ncol,this%nlev  ))
    allocate(this%lnp_lev   (this%ncol,this%nlev+1))
    allocate(this%omg       (this%ncol,this%nlev  ))
    ! Height
    allocate(this%z         (this%ncol,this%nlev  ))
    allocate(this%z_lev     (this%ncol,this%nlev+1))
    allocate(this%dz        (this%ncol,this%nlev  ))
    ! Tracers
  if (ntracers > 0) then
    allocate(this%q(this%ncol,this%nlev,ntracers))
  end if
    ! Moisture
    if (idx_qv /= 0) this%qv => this%q(:,:,idx_qv)
    if (idx_qc /= 0) this%qc => this%q(:,:,idx_qc)
    if (idx_qi /= 0) this%qi => this%q(:,:,idx_qi)
    if (idx_qr /= 0) this%qr => this%q(:,:,idx_qr)
    if (idx_qs /= 0) this%qs => this%q(:,:,idx_qs)
    if (idx_qg /= 0) this%qg => this%q(:,:,idx_qg)
    if (idx_qh /= 0) this%qh => this%q(:,:,idx_qh)
    ! Ozone
    if (idx_qo3 /= 0) this%qo3 => this%q(:,:,idx_qo3)
    ! Stability
    allocate(this%n2_lev    (this%ncol,this%nlev+1))
    allocate(this%ri_lev    (this%ncol,this%nlev+1))
    ! Surface layer
    allocate(this%emis      (this%ncol            ))
    allocate(this%alb       (this%ncol            ))
    allocate(this%ps        (this%ncol            ))
    allocate(this%ts        (this%ncol            ))
    allocate(this%land      (this%ncol            ))
    allocate(this%hfx       (this%ncol            ))
    allocate(this%qfx       (this%ncol            ))
    allocate(this%z0        (this%ncol            ))
    allocate(this%ustar     (this%ncol            ))
    allocate(this%ptstar    (this%ncol            ))
    allocate(this%psim      (this%ncol            ))
    allocate(this%psih      (this%ncol            ))
    allocate(this%rib       (this%ncol            ))
    allocate(this%uos       (this%ncol            ))
    allocate(this%vos       (this%ncol            ))
    ! Boundary layer
    allocate(this%pblh      (this%ncol            ))
    allocate(this%pblk      (this%ncol            ))
    allocate(this%wstar     (this%ncol            ))
    allocate(this%exch_h    (this%ncol,this%nlev  ))
    allocate(this%delta     (this%ncol,this%nlev  ))
    ! Precipitation
    allocate(this%precl     (this%ncol            ))
    ! Mars
    select case (planet)
    case ('mars')
      allocate(this%co2ice  (this%ncol            ))
    end select
    ! Others
    allocate(this%rho       (this%ncol,this%nlev  ))
    allocate(this%cp        (this%ncol,this%nlev  ))
    allocate(this%cv        (this%ncol,this%nlev  ))
    allocate(this%tep       (this%ncol,this%nlev  ))
    allocate(this%tev       (this%ncol,this%nlev  ))

    icol = 0
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        icol = icol + 1
        this%i   (icol) = i
        this%j   (icol) = j
        this%lon (icol) = mesh%full_lon_deg(i)
        this%lat (icol) = mesh%full_lat_deg(j)
        this%area(icol) = mesh%area_cell(j)
      end do
    end do

  end subroutine pstate_init

  subroutine pstate_clear(this)

    class(pstate_type), intent(inout) :: this

    if (allocated(this%i        )) deallocate(this%i        )
    if (allocated(this%j        )) deallocate(this%j        )
    if (allocated(this%lon      )) deallocate(this%lon      )
    if (allocated(this%lat      )) deallocate(this%lat      )
    ! Wind
    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%wsb      )) deallocate(this%wsb      )
    if (allocated(this%u10      )) deallocate(this%u10      )
    if (allocated(this%v10      )) deallocate(this%v10      )
    ! Temperature
    if (allocated(this%t        )) deallocate(this%t        )
    if (allocated(this%tv       )) deallocate(this%tv       )
    if (allocated(this%pt       )) deallocate(this%pt       )
    if (allocated(this%ptv      )) deallocate(this%ptv      )
    if (allocated(this%t_sfc    )) deallocate(this%t_sfc    )
    ! Pressure
    if (allocated(this%p        )) deallocate(this%p        )
    if (allocated(this%p_lev    )) deallocate(this%p_lev    )
    if (allocated(this%pk       )) deallocate(this%pk       )
    if (allocated(this%pk_lev   )) deallocate(this%pk_lev   )
    if (allocated(this%dp       )) deallocate(this%dp       )
    if (allocated(this%rdp      )) deallocate(this%rdp      )
    if (allocated(this%lnp_lev  )) deallocate(this%lnp_lev  )
    ! Height
    if (allocated(this%z        )) deallocate(this%z        )
    if (allocated(this%z_lev    )) deallocate(this%z_lev    )
    if (allocated(this%dz       )) deallocate(this%dz       )
    ! Tracers
    if (allocated(this%q        )) deallocate(this%q        )
    ! Stability
    if (allocated(this%n2_lev   )) deallocate(this%n2_lev   )
    if (allocated(this%ri_lev   )) deallocate(this%ri_lev   )
    ! Surface layer
    if (allocated(this%emis     )) deallocate(this%emis     )
    if (allocated(this%alb      )) deallocate(this%alb      )
    if (allocated(this%ps       )) deallocate(this%ps       )
    if (allocated(this%ts       )) deallocate(this%ts       )
    if (allocated(this%land     )) deallocate(this%land     )
    if (allocated(this%hfx      )) deallocate(this%hfx      )
    if (allocated(this%qfx      )) deallocate(this%qfx      )
    if (allocated(this%z0       )) deallocate(this%z0       )
    if (allocated(this%ustar    )) deallocate(this%ustar    )
    if (allocated(this%ptstar   )) deallocate(this%ptstar   )
    if (allocated(this%psim     )) deallocate(this%psim     )
    if (allocated(this%psih     )) deallocate(this%psih     )
    if (allocated(this%rib      )) deallocate(this%rib      )
    if (allocated(this%uos      )) deallocate(this%uos      )
    if (allocated(this%vos      )) deallocate(this%vos      )
    ! Boundary layer
    if (allocated(this%pblh     )) deallocate(this%pblh     )
    if (allocated(this%pblk     )) deallocate(this%pblk     )
    if (allocated(this%wstar    )) deallocate(this%wstar    )
    if (allocated(this%exch_h   )) deallocate(this%exch_h   )
    if (allocated(this%delta    )) deallocate(this%delta    )
    ! Precipitation
    if (allocated(this%precl    )) deallocate(this%precl    )
    ! Mars
    if (allocated(this%co2ice   )) deallocate(this%co2ice   )
    ! Others
    if (allocated(this%rho      )) deallocate(this%rho      )
    if (allocated(this%cp       )) deallocate(this%cp       )
    if (allocated(this%cv       )) deallocate(this%cv       )
    if (allocated(this%tep      )) deallocate(this%tep      )
    if (allocated(this%tev      )) deallocate(this%tev      )

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
    allocate(this%dqdt      (this%ncol,this%nlev,ntracers))
    allocate(this%dptdt     (this%ncol,this%nlev))
    allocate(this%dptdt_rad (this%ncol,this%nlev)); this%dptdt_rad = 0

  end subroutine ptend_init

  subroutine ptend_clear(this)

    class(ptend_type), intent(inout) :: this

    if (allocated(this%dudt     )) deallocate(this%dudt     )
    if (allocated(this%dvdt     )) deallocate(this%dvdt     )
    if (allocated(this%dtdt     )) deallocate(this%dtdt     )
    if (allocated(this%dqdt     )) deallocate(this%dqdt     )
    if (allocated(this%dptdt    )) deallocate(this%dptdt    )
    if (allocated(this%dptdt_rad)) deallocate(this%dptdt_rad)

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
    this%dqdt  = 0
    this%updated_qv = .false.
    this%updated_qc = .false.
    this%updated_qi = .false.

  end subroutine ptend_reset

end module physics_types_mod
