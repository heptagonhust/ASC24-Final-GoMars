module dynamics_types_mod

  use const_mod
  use namelist_mod
  use latlon_mesh_mod
  use allocator_mod
  use tracer_types_mod

  implicit none

  private

  public dstate_type
  public dtend_type
  public static_type
  public aux_array_type

  ! NOTE:
  !   Variables with '_lon', '_lat' and '_lev' are on the half grids on the corresponding direction,
  !   and '_p' indicates that the variable is perturbed.
  type dstate_type
    type(latlon_mesh_type), pointer :: mesh => null()
    ! For nesting
    integer :: id = 0
    type(dstate_type), pointer :: parent => null()
    real(r8), allocatable, dimension(:,:,:) :: u                 ! Zonal wind speed at cell center (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: v                 ! Meridional wind speed at cell center (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: u_lon             ! Zonal wind speed at lon edge (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: v_lat             ! Meridional wind speed at lat edge (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: we_lev            ! Vertical coordinate speed multiplied by 𝛛π/𝛛η
    real(r8), allocatable, dimension(:,:,:) :: gz                ! Geopotential (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: gz_lev            ! Geopotential height on half levels (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: dmg               ! Mass
    real(r8), allocatable, dimension(:,:,:) :: dmg_lev           ! Mass on half levels
    real(r8), allocatable, dimension(:,:,:) :: pt                ! Potential temperature
    real(r8), allocatable, dimension(:,:,:) :: t                 ! Temperature
    real(r8), allocatable, dimension(:,:,:) :: tv                ! Virtual temperature
    real(r8), allocatable, dimension(:,:,:) :: mg                ! Dry air weight on full levels
    real(r8), allocatable, dimension(:,:,:) :: mg_lev            ! Dry air weight on half levels
    real(r8), allocatable, dimension(:,:  ) :: mgs               ! Surface dry air weight
    real(r8), allocatable, dimension(:,:,:) :: ph                ! Hydrostatic pressure (dry air and water vapor) on full levels
    real(r8), allocatable, dimension(:,:,:) :: ph_lev            ! Hydrostatic pressure (dry air and water vapor) on half levels
    real(r8), pointer    , dimension(:,:  ) :: phs               ! Surface hydrostatic pressure (dry air and water vapor)
    real(r8), allocatable, dimension(:,:,:) :: rhod              ! Dry air density
    ! Nonhydrostatic variables
    real(r8), allocatable, dimension(:,:,:) :: we
    real(r8), allocatable, dimension(:,:,:) :: w                 ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev             ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev_lon         ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev_lat         ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: gz_lev_lon        ! Geopotential
    real(r8), allocatable, dimension(:,:,:) :: gz_lev_lat        ! Geopotential
    real(r8), pointer    , dimension(:,:,:) :: p                 ! Full pressure on full levels
    real(r8), pointer    , dimension(:,:,:) :: p_lev             ! Full pressure on half levels
    real(r8), allocatable, dimension(:,:  ) :: ps                ! Surface full pressure
    real(r8), allocatable, dimension(:,:,:) :: p_lev_lon         ! Full pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: p_lev_lat         ! Full pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: u_lev_lon
    real(r8), allocatable, dimension(:,:,:) :: v_lev_lat
    real(r8), allocatable, dimension(:,:,:) :: mf_lev_lon_n      ! Mass flux on zonal edge and half level
    real(r8), allocatable, dimension(:,:,:) :: mf_lev_lat_n      ! Mass flux on merdional edge and half level
    ! Total diagnostics
    real(r8) tm
    real(r8) te, te_ke, te_ie, te_pe
    real(r8) tpe
  contains
    procedure :: init         => dstate_init
    procedure :: clear        => dstate_clear
    procedure :: c2a          => dstate_c2a
    procedure :: a2c          => dstate_a2c
    generic :: operator(+)    => dstate_add
    generic :: operator(*)    => dstate_mul
    generic :: operator(/)    => dstate_div
    generic :: assignment(=)  => dstate_assign
    procedure, pass(x) :: dstate_add, dstate_mul, dstate_div, dstate_assign
    final :: dstate_final
  end type dstate_type

  type dtend_type
    type(latlon_mesh_type), pointer :: filter_mesh => null()
    type(latlon_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:,:) :: du
    real(r8), allocatable, dimension(:,:,:) :: dv
    real(r8), allocatable, dimension(:,:,:) :: dgz
    real(r8), allocatable, dimension(:,:,:) :: dpt
    real(r8), allocatable, dimension(:,:  ) :: dmgs
#ifdef OUTPUT_H1_DTEND
    real(r8), allocatable, dimension(:,:,:) :: dudt_coriolis
    real(r8), allocatable, dimension(:,:,:) :: dvdt_coriolis
    real(r8), allocatable, dimension(:,:,:) :: dudt_wedudeta
    real(r8), allocatable, dimension(:,:,:) :: dvdt_wedvdeta
    real(r8), allocatable, dimension(:,:,:) :: dudt_dkedx
    real(r8), allocatable, dimension(:,:,:) :: dvdt_dkedy
    real(r8), allocatable, dimension(:,:,:) :: dudt_pgf
    real(r8), allocatable, dimension(:,:,:) :: dvdt_pgf
#endif
    logical :: update_u   = .false.
    logical :: update_v   = .false.
    logical :: update_gz  = .false.
    logical :: update_pt  = .false.
    logical :: update_mgs = .false.
    logical :: copy_gz    = .false.
    logical :: copy_pt    = .false.
    logical :: copy_mgs   = .false.
    ! Nonhydrostatic tendencies
    real(r8), allocatable, dimension(:,:,:) :: adv_gz
    real(r8), allocatable, dimension(:,:,:) :: adv_w
  contains
    procedure :: init         => dtend_init
    procedure :: reset_flags  => dtend_reset_flags
    procedure :: clear        => dtend_clear
    generic :: operator(+)    => dtend_add
    generic :: operator(*)    => dtend_mul
    generic :: operator(/)    => dtend_div
    generic :: assignment(=)  => dtend_assign
    procedure, pass(x) :: dtend_add, dtend_mul, dtend_div, dtend_assign
    final :: dtend_final
  end type dtend_type

  type static_type
    type(latlon_mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:) :: landmask
    ! Topography
    real(r8), allocatable, dimension(:,:) :: gzs
    real(r8), allocatable, dimension(:,:) :: zs_std
    real(r8), allocatable, dimension(:,:) :: dzsdlon
    real(r8), allocatable, dimension(:,:) :: dzsdlat
    ! Reference surface pressure
    real(r8), allocatable, dimension(:,:) :: ref_ps
    real(r8), allocatable, dimension(:,:) :: ref_ps_smth
    real(r8), allocatable, dimension(:,:) :: ref_ps_perb
    ! Coriolis parameters
    real(8), allocatable, dimension(:  ) :: full_f
    real(8), allocatable, dimension(:  ) :: half_f
    ! Weight for constructing tangential wind
    real(8), allocatable, dimension(:,:) :: full_tangent_wgt
    real(8), allocatable, dimension(:,:) :: half_tangent_wgt
  contains
    procedure :: init_stage1 => static_init_stage1
    procedure :: init_stage2 => static_init_stage2
    procedure :: clear       => static_clear
    final :: static_final
  end type static_type

  type aux_array_type
    type(latlon_mesh_type), pointer :: filter_mesh => null()
    type(latlon_mesh_type), pointer :: mesh => null()
    ! Smagorinsky damping variables
    real(r8), allocatable, dimension(:,:,:) :: smag_t            ! tension strain
    real(r8), allocatable, dimension(:,:,:) :: smag_s            ! shear strain on vertex
    real(r8), allocatable, dimension(:,:,:) :: kmh               ! nonlinear diffusion coef
    real(r8), allocatable, dimension(:,:,:) :: kmh_lon           ! nonlinear diffusion coef on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: kmh_lat           ! nonlinear diffusion coef on meridional edge
    real(r8), allocatable, dimension(:,:,:) :: v_lon             ! Meridional wind speed at lon edge (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: u_lat             ! Zonal wind speed at lat edge (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: ke                ! Kinetic energy
    real(r8), allocatable, dimension(:,:,:) :: pv_lon            ! Potential vorticity on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: pv_lat            ! Potential vorticity on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: dmg_lon           ! Mass on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: dmg_lat           ! Mass on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: dmg_vtx           ! Mass on vertex
    real(r8), allocatable, dimension(:,:,:) :: pkh_lev           ! Exner pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: we_lev_lon        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: we_lev_lat        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: ptf_lon           ! Potential temperature on the zonal edge
    real(r8), allocatable, dimension(:,:,:) :: ptf_lat           ! Potential temperature on the merdional edge
    real(r8), allocatable, dimension(:,:,:) :: ptf_lev           ! Potential temperature on the vertical edge
    real(r8), allocatable, dimension(:,:,:) :: mfx_lon           ! Normal mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mfy_lat           ! Normal mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: mfx_lat           ! Tangient mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mfy_lon           ! Tangient mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: vor               ! Vorticity (s-1)
    real(r8), allocatable, dimension(:,:,:) :: pv                ! Potential vorticity
    real(r8), allocatable, dimension(:,:,:) :: div               ! Divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: div2              ! Laplacian of divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: dmf
    real(r8), allocatable, dimension(:,:,:) :: omg               ! Vertical pressure velocity (Pa s-1)
    ! Tendencies from Smagorinsky diffusion
    real(r8), allocatable, dimension(:,:,:  ) :: dudt_damp
    real(r8), allocatable, dimension(:,:,:  ) :: dvdt_damp
    real(r8), allocatable, dimension(:,:,:  ) :: dptdt_damp
    real(r8), allocatable, dimension(:,:,:,:) :: dqdt_damp
    ! Tendencies from physics
    real(r8), allocatable, dimension(:,:,:  ) :: dudt_phys
    real(r8), allocatable, dimension(:,:,:  ) :: dvdt_phys
    real(r8), allocatable, dimension(:,:,:  ) :: dptdt_phys
    real(r8), allocatable, dimension(:,:,:,:) :: dqdt_phys
    ! Perturbed quantities for calculating HPGF
    real(r8), allocatable, dimension(:,:,:  ) :: p_ptb
    real(r8), allocatable, dimension(:,:,:  ) :: gz_ptb
    real(r8), allocatable, dimension(:,:,:  ) :: dp_ptb
    real(r8), allocatable, dimension(:,:,:  ) :: ad_ptb
  contains
    procedure :: init      => aux_array_init
    procedure :: init_phys => aux_array_init_phys
    procedure :: clear     => aux_array_clear
    final :: aux_array_final
  end type aux_array_type

contains

  subroutine dstate_init(this, filter_mesh, mesh)

    class(dstate_type), intent(inout), target :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u                , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v                , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u_lon            , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v_lat            , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%we_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%gz               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%dmg              , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmg_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%t                , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%tv               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mg               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mg_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%mgs              , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%ph               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ph_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%rhod             , full_lon=.true., full_lat=.true., full_lev=.true.)

    call allocate_array(filter_mesh, this%pt        , full_lon=.true., full_lat=.true., full_lev=.true.)

    if (baroclinic) then
      this%phs(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme) => this%ph_lev(:,:,mesh%half_kde)
    end if

    if (nonhydrostatic) then
      call allocate_array(mesh, this%we             , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%w              , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%w_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%w_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%w_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%gz_lev_lon     , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%gz_lev_lat     , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%p              , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%p_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%p_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%p_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%u_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%v_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%mf_lev_lon_n   , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%mf_lev_lat_n   , full_lon=.true., half_lat=.true., half_lev=.true.)
    else
      this%p     => this%ph
      this%p_lev => this%ph_lev
    end if
    call allocate_array(mesh, this%ps, full_lon=.true., full_lat=.true.)

  end subroutine dstate_init

  subroutine dstate_clear(this)

    class(dstate_type), intent(inout) :: this

    if (allocated(this%u                )) deallocate(this%u                )
    if (allocated(this%v                )) deallocate(this%v                )
    if (allocated(this%u_lon            )) deallocate(this%u_lon            )
    if (allocated(this%v_lat            )) deallocate(this%v_lat            )
    if (allocated(this%we_lev           )) deallocate(this%we_lev           )
    if (allocated(this%gz               )) deallocate(this%gz               )
    if (allocated(this%gz_lev           )) deallocate(this%gz_lev           )
    if (allocated(this%dmg              )) deallocate(this%dmg              )
    if (allocated(this%dmg_lev          )) deallocate(this%dmg_lev          )
    if (allocated(this%pt               )) deallocate(this%pt               )
    if (allocated(this%t                )) deallocate(this%t                )
    if (allocated(this%tv               )) deallocate(this%tv               )
    if (allocated(this%mg               )) deallocate(this%mg               )
    if (allocated(this%mg_lev           )) deallocate(this%mg_lev           )
    if (allocated(this%mgs              )) deallocate(this%mgs              )
    if (allocated(this%ph               )) deallocate(this%ph               )
    if (allocated(this%ph_lev           )) deallocate(this%ph_lev           )
    if (allocated(this%rhod             )) deallocate(this%rhod             )

    if (allocated(this%we               )) deallocate(this%we               )
    if (allocated(this%w                )) deallocate(this%w                )
    if (allocated(this%w_lev            )) deallocate(this%w_lev            )
    if (allocated(this%w_lev_lon        )) deallocate(this%w_lev_lon        )
    if (allocated(this%w_lev_lat        )) deallocate(this%w_lev_lat        )
    if (allocated(this%gz_lev_lon       )) deallocate(this%gz_lev_lon       )
    if (allocated(this%gz_lev_lat       )) deallocate(this%gz_lev_lat       )
    if (allocated(this%p_lev_lon        )) deallocate(this%p_lev_lon        )
    if (allocated(this%p_lev_lat        )) deallocate(this%p_lev_lat        )
    if (allocated(this%u_lev_lon        )) deallocate(this%u_lev_lon        )
    if (allocated(this%v_lev_lat        )) deallocate(this%v_lev_lat        )
    if (allocated(this%mf_lev_lon_n     )) deallocate(this%mf_lev_lon_n     )
    if (allocated(this%mf_lev_lat_n     )) deallocate(this%mf_lev_lat_n     )

    if (nonhydrostatic) then
      if (associated(this%p             )) deallocate(this%p                )
      if (associated(this%p_lev         )) deallocate(this%p_lev            )
    end if

  end subroutine dstate_clear

  subroutine dstate_final(this)

    type(dstate_type), intent(inout) :: this

    call this%clear()

  end subroutine dstate_final

  subroutine dstate_a2c(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    do k = this%mesh%full_kds, this%mesh%full_kde
      do j = this%mesh%full_jds_no_pole, this%mesh%full_jde_no_pole
        do i = this%mesh%half_ids, this%mesh%half_ide
          this%u_lon(i,j,k) = 0.5_r8 * (this%u(i,j,k) + this%u(i+1,j,k))
        end do
      end do
      do j = this%mesh%half_jds, this%mesh%half_jde
        do i = this%mesh%full_ids, this%mesh%full_ide
          this%v_lat(i,j,k) = 0.5_r8 * (this%v(i,j,k) + this%v(i,j+1,k))
        end do
      end do
    end do

  end subroutine dstate_a2c

  subroutine dstate_c2a(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    do k = this%mesh%full_kds, this%mesh%full_kde
      do j = this%mesh%full_jds_no_pole, this%mesh%full_jde_no_pole
        do i = this%mesh%full_ids, this%mesh%full_ide
          this%u(i,j,k) = 0.5_r8 * (this%u_lon(i,j,k) + this%u_lon(i-1,j,k))
          this%v(i,j,k) = 0.5_r8 * (this%v_lat(i,j,k) + this%v_lat(i,j-1,k))
        end do
      end do
    end do

  end subroutine dstate_c2a

  function dstate_add(x, y) result(res)

    class(dstate_type), intent(in) :: x
    class(dstate_type), intent(in) :: y

    type(dstate_type) res

    if (hydrostatic) then
      res%u_lon = x%u_lon + y%u_lon
      res%v_lat = x%v_lat + y%v_lat
      res%pt    = x%pt    + y%pt
      res%mgs   = x%mgs   + y%mgs
    else if (nonhydrostatic) then
    else
      res%u_lon = x%u_lon + y%u_lon
      res%v_lat = x%v_lat + y%v_lat
      res%gz    = x%gz    + y%gz
    end if

  end function dstate_add

  function dstate_mul(s, x) result(res)

    real(r8), intent(in) :: s
    class(dstate_type), intent(in) :: x

    type(dstate_type) res

    if (hydrostatic) then
      res%u_lon = s * x%u_lon
      res%v_lat = s * x%v_lat
      res%pt    = s * x%pt
      res%mgs   = s * x%mgs
    else if (nonhydrostatic) then
    else
      res%u_lon = s * x%u_lon
      res%v_lat = s * x%v_lat
      res%gz    = s * x%gz
    end if

  end function dstate_mul

  function dstate_div(x, s) result(res)

    class(dstate_type), intent(in) :: x
    real(r8), intent(in) :: s

    type(dstate_type) res

    if (hydrostatic) then
      res%u_lon = x%u_lon / s
      res%v_lat = x%v_lat / s
      res%pt    = x%pt / s
      res%mgs   = x%mgs / s
    else if (nonhydrostatic) then
    else
      res%u_lon = x%u_lon / s
      res%v_lat = x%v_lat / s
      res%gz    = x%gz / s
    end if

  end function dstate_div

  subroutine dstate_assign(x, y)

    class(dstate_type), intent(inout) :: x
    class(dstate_type), intent(in) :: y

    if (hydrostatic) then
      x%u_lon = y%u_lon
      x%v_lat = y%v_lat
      x%pt    = y%pt
      x%mgs   = y%mgs
    else if (nonhydrostatic) then
    else
      x%u_lon = y%u_lon
      x%v_lat = y%v_lat
      x%gz    = y%gz
    end if

  end subroutine dstate_assign

  subroutine dtend_init(this, filter_mesh, mesh)

    class(dtend_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_mesh_type), intent(in), target :: mesh

    call this%clear()

    this%filter_mesh => filter_mesh
    this%mesh => mesh

    call allocate_array(filter_mesh, this%du  , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(filter_mesh, this%dv  , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(filter_mesh, this%dpt , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(filter_mesh, this%dmgs, full_lon=.true., full_lat=.true.                 )

    if (nonhydrostatic .or. .not. baroclinic) then
      call allocate_array(filter_mesh, this%dgz , full_lon=.true., full_lat=.true., full_lev=.true.)
    end if
    if (nonhydrostatic) then
      call allocate_array(mesh, this%adv_gz, full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%adv_w , full_lon=.true., full_lat=.true., half_lev=.true.)
    end if

#ifdef OUTPUT_H1_DTEND
    call allocate_array(mesh, this%dudt_coriolis, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvdt_coriolis, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dudt_wedudeta, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvdt_wedvdeta, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dudt_dkedx   , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvdt_dkedy   , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dudt_pgf     , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvdt_pgf     , full_lon=.true., half_lat=.true., full_lev=.true.)
#endif

  end subroutine dtend_init

  subroutine dtend_reset_flags(this)

    class(dtend_type), intent(inout) :: this

    this%du = 0
    this%dv = 0

    this%update_u   = .false.
    this%update_v   = .false.
    this%update_gz  = .false.
    this%update_pt  = .false.
    this%update_mgs = .false.
    this%copy_gz    = .false.
    this%copy_pt    = .false.
    this%copy_mgs   = .false.

  end subroutine dtend_reset_flags

  subroutine dtend_clear(this)

    class(dtend_type), intent(inout) :: this

    if (allocated(this%du      )) deallocate(this%du      )
    if (allocated(this%dv      )) deallocate(this%dv      )
    if (allocated(this%dgz     )) deallocate(this%dgz     )
    if (allocated(this%dpt     )) deallocate(this%dpt     )
    if (allocated(this%dmgs    )) deallocate(this%dmgs    )

    if (allocated(this%adv_gz  )) deallocate(this%adv_gz  )
    if (allocated(this%adv_w   )) deallocate(this%adv_w   )

#ifdef OUTPUT_H1_DTEND
    if (allocated(this%dudt_coriolis)) deallocate(this%dudt_coriolis)
    if (allocated(this%dvdt_coriolis)) deallocate(this%dvdt_coriolis)
    if (allocated(this%dudt_wedudeta)) deallocate(this%dudt_wedudeta)
    if (allocated(this%dvdt_wedvdeta)) deallocate(this%dvdt_wedvdeta)
    if (allocated(this%dudt_dkedx   )) deallocate(this%dudt_dkedx   )
    if (allocated(this%dvdt_dkedy   )) deallocate(this%dvdt_dkedy   )
    if (allocated(this%dudt_pgf     )) deallocate(this%dudt_pgf     )
    if (allocated(this%dvdt_pgf     )) deallocate(this%dvdt_pgf     )
#endif

  end subroutine dtend_clear

  subroutine dtend_final(this)

    type(dtend_type), intent(inout) :: this

    call this%clear()

  end subroutine dtend_final

  function dtend_add(x, y) result(res)

    class(dtend_type), intent(in) :: x
    class(dtend_type), intent(in) :: y

    type(dtend_type) res

    if (x%update_u .and. y%update_u) then
      res%du = x%du + y%du
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v .and. y%update_v) then
      res%dv = x%dv + y%dv
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_mgs .and. y%update_mgs) then
        res%dmgs = x%dmgs + y%dmgs
        res%update_mgs = .true.
      else
        res%update_mgs = .false.
      end if
      if (x%update_pt .and. y%update_pt) then
        res%dpt = x%dpt + y%dpt
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz .and. y%update_gz) then
      res%dgz = x%dgz + y%dgz
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function dtend_add

  function dtend_mul(s, x) result(res)

    real(r8), intent(in) :: s
    class(dtend_type), intent(in) :: x

    type(dtend_type) res

    if (x%update_u) then
      res%du = s * x%du
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v) then
      res%dv = s * x%dv
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_mgs) then
        res%dmgs = s * x%dmgs
        res%update_mgs = .true.
      else
        res%update_mgs = .false.
      end if
      if (x%update_pt) then
        res%dpt = s * x%dpt
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz) then
      res%dgz = s * x%dgz
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function dtend_mul

  function dtend_div(x, s) result(res)

    class(dtend_type), intent(in) :: x
    real(r8), intent(in) :: s

    type(dtend_type) res

    if (x%update_u) then
      res%du = x%du / s
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v) then
      res%dv = x%dv / s
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_mgs) then
        res%dmgs = x%dmgs / s
        res%update_mgs = .true.
      else
        res%update_mgs = .false.
      end if
      if (x%update_pt) then
        res%dpt = x%dpt / s
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz) then
      res%dgz = x%dgz / s
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function dtend_div

  subroutine dtend_assign(x, y)

    class(dtend_type), intent(inout) :: x
    class(dtend_type), intent(in) :: y

    if (y%update_u) then
      x%du = y%du
      x%update_u = .true.
    else
      x%update_u = .false.
    end if
    if (y%update_v) then
      x%dv = y%dv
      x%update_v = .true.
    else
      x%update_v = .false.
    end if
    if (baroclinic) then
      if (y%update_mgs) then
        x%dmgs = y%dmgs
        x%update_mgs = .true.
      else
        x%update_mgs = .false.
      end if
      if (y%update_pt) then
        x%dpt = y%dpt
        x%update_pt = .true.
      else
        x%update_pt = .false.
      end if
    else if (y%update_gz) then
      x%dgz = y%dgz
      x%update_gz = .true.
    else
      x%update_gz = .false.
    end if

  end subroutine dtend_assign

  subroutine static_init_stage1(this, filter_mesh, mesh)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(filter_mesh, this%gzs, full_lon=.true., full_lat=.true.)

    call allocate_array(mesh, this%landmask   , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%zs_std     , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlon    , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlat    , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ref_ps     , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%ref_ps_smth, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%ref_ps_perb, full_lon=.true., full_lat=.true.)

    allocate(this%full_f             (this%mesh%full_jms:this%mesh%full_jme)); this%full_f              = inf
    allocate(this%half_f             (this%mesh%half_jms:this%mesh%half_jme)); this%half_f              = inf
    allocate(this%full_tangent_wgt (2,this%mesh%full_jms:this%mesh%full_jme)); this%full_tangent_wgt    = inf
    allocate(this%half_tangent_wgt (2,this%mesh%half_jms:this%mesh%half_jme)); this%half_tangent_wgt    = inf

  end subroutine static_init_stage1

  subroutine static_init_stage2(this)

    class(static_type), intent(inout) :: this

    integer j

    do j = this%mesh%full_jds, this%mesh%full_jde
      this%full_f(j) = 2 * omega * this%mesh%full_sin_lat(j)
    end do
    do j = this%mesh%half_jds, this%mesh%half_jde
      this%half_f(j) = 2 * omega * this%mesh%half_sin_lat(j)
    end do


    !  ____________________                 ____________________                  ____________________                  ____________________
    ! |          |         |               |          |         |                |          |         |                |          |         |
    ! |          |         |               |          |         |                |          |         |                |          |         |
    ! |          |         |               |          |         |                |          |         |                |          |         |
    ! |          |         |               |          |         |                |          |         |                |          |         |
    ! |_____o____|____o____|   j           |_____o____|____*____|   j            |_____*____|____o____|   j            |_____o____|____o____|   j
    ! |          |////|////|               |          |////|    |                |/////|    |         |                |     |    |         |
    ! |          |/3//|/2//|               |          |////|    |                |/////|    |         |                |     |    |         |
    ! |          x---------|   j           |          x---------|   j            |-----|----x         |   j            |-----|----x         |   j
    ! |          |    |/1//|               |          |    |    |                |/////|////|         |                |     |////|         |
    ! |_____o____|____*____|   j - 1       |_____o____|____o____|   j - 1        |_____o____|____o____|   j - 1        |_____*____|____o____|   j - 1
    !       i    i   i+1                         i    i   i+1                          i    i   i+1                          i
    !
    !
    !       [ 1    As_1 + As_2 + As_3]
    ! w = - [--- - ------------------]
    !       [ 2        A_{i+1,j}     ]
    !
    !

    select case (tangent_wgt_scheme)
    case ('classic')
      do j = this%mesh%full_jds_no_pole, this%mesh%full_jde_no_pole
        this%full_tangent_wgt(1,j) = this%mesh%le_lat(j-1) / this%mesh%de_lon(j) * 0.25d0
        this%full_tangent_wgt(2,j) = this%mesh%le_lat(j  ) / this%mesh%de_lon(j) * 0.25d0
      end do

      do j = this%mesh%half_jds, this%mesh%half_jde
        this%half_tangent_wgt(1,j) = this%mesh%le_lon(j  ) / this%mesh%de_lat(j) * 0.25d0
        this%half_tangent_wgt(2,j) = this%mesh%le_lon(j+1) / this%mesh%de_lat(j) * 0.25d0
      end do
    case ('th09')
      do j = this%mesh%full_jds_no_pole, this%mesh%full_jde_no_pole
        this%full_tangent_wgt(1,j) = this%mesh%le_lat(j-1) / this%mesh%de_lon(j) * this%mesh%area_subcell(2,j  ) / this%mesh%area_cell(j  )
        this%full_tangent_wgt(2,j) = this%mesh%le_lat(j  ) / this%mesh%de_lon(j) * this%mesh%area_subcell(1,j  ) / this%mesh%area_cell(j  )
      end do

      do j = this%mesh%half_jds, this%mesh%half_jde
        this%half_tangent_wgt(1,j) = this%mesh%le_lon(j  ) / this%mesh%de_lat(j) * this%mesh%area_subcell(1,j  ) / this%mesh%area_cell(j  )
        this%half_tangent_wgt(2,j) = this%mesh%le_lon(j+1) / this%mesh%de_lat(j) * this%mesh%area_subcell(2,j+1) / this%mesh%area_cell(j+1)
      end do
    end select

  end subroutine static_init_stage2

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%landmask        )) deallocate(this%landmask        )
    if (allocated(this%gzs             )) deallocate(this%gzs             )
    if (allocated(this%zs_std          )) deallocate(this%zs_std          )
    if (allocated(this%dzsdlon         )) deallocate(this%dzsdlon         )
    if (allocated(this%dzsdlat         )) deallocate(this%dzsdlat         )
    if (allocated(this%ref_ps          )) deallocate(this%ref_ps          )
    if (allocated(this%ref_ps_smth     )) deallocate(this%ref_ps_smth     )
    if (allocated(this%ref_ps_perb     )) deallocate(this%ref_ps_perb     )
    if (allocated(this%full_f          )) deallocate(this%full_f          )
    if (allocated(this%half_f          )) deallocate(this%half_f          )
    if (allocated(this%full_tangent_wgt)) deallocate(this%full_tangent_wgt)
    if (allocated(this%half_tangent_wgt)) deallocate(this%half_tangent_wgt)

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

  subroutine aux_array_init(this, filter_mesh, mesh)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_mesh_type), intent(in), target :: mesh

    this%filter_mesh => filter_mesh
    this%mesh => mesh

    if (use_smag_damp) then
      call allocate_array(mesh, this%smag_t       , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%smag_s       , half_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh          , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh_lon      , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh_lat      , full_lon=.true., half_lat=.true., full_lev=.true.)
    end if
    call allocate_array(mesh, this%v_lon          , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u_lat          , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ke             , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lon         , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lat         , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmg_lon        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmg_lat        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmg_vtx        , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pkh_lev        , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%we_lev_lon     , half_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%we_lev_lat     , full_lon=.true., half_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%ptf_lon        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ptf_lat        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ptf_lev        , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%mfx_lon        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfy_lat        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfx_lat        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfy_lon        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%vor            , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv             , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%div            , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%div2           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmf            , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%omg            , full_lon=.true., full_lat=.true., full_lev=.true.)

    call allocate_array(this%mesh, this% dudt_damp, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this% dvdt_damp, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%dptdt_damp, full_lon=.true., full_lat=.true., full_lev=.true.)

    if (pgf_scheme == 'ptb') then
      call allocate_array(mesh, this%p_ptb        , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%gz_ptb       , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%dp_ptb       , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%ad_ptb       , full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine aux_array_init

  subroutine aux_array_init_phys(this)

    class(aux_array_type), intent(inout) :: this

    call allocate_array(this%mesh, this% dqdt_damp, full_lon=.true., full_lat=.true., full_lev=.true., extra_dim=ntracers)

    if (trim(physics_suite) /= '') then
      if (filter_ptend) then
        call allocate_array(this%filter_mesh, this% dudt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%filter_mesh, this% dvdt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%filter_mesh, this%dptdt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%filter_mesh, this% dqdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., extra_dim=ntracers)
      else
        call allocate_array(this%mesh, this% dudt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this% dvdt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%dptdt_phys, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this% dqdt_phys, full_lon=.true., full_lat=.true., full_lev=.true., extra_dim=ntracers)
      end if
    end if

  end subroutine aux_array_init_phys

  subroutine aux_array_clear(this)

    class(aux_array_type), intent(inout) :: this

    if (allocated(this%smag_t           )) deallocate(this%smag_t           )
    if (allocated(this%smag_s           )) deallocate(this%smag_s           )
    if (allocated(this%kmh              )) deallocate(this%kmh              )
    if (allocated(this%kmh_lon          )) deallocate(this%kmh_lon          )
    if (allocated(this%kmh_lat          )) deallocate(this%kmh_lat          )
    if (allocated(this%v_lon            )) deallocate(this%v_lon            )
    if (allocated(this%u_lat            )) deallocate(this%u_lat            )
    if (allocated(this%ke               )) deallocate(this%ke               )
    if (allocated(this%pv_lon           )) deallocate(this%pv_lon           )
    if (allocated(this%pv_lat           )) deallocate(this%pv_lat           )
    if (allocated(this%dmg_lon          )) deallocate(this%dmg_lon          )
    if (allocated(this%dmg_lat          )) deallocate(this%dmg_lat          )
    if (allocated(this%dmg_vtx          )) deallocate(this%dmg_vtx          )
    if (allocated(this%pkh_lev          )) deallocate(this%pkh_lev          )
    if (allocated(this%we_lev_lon       )) deallocate(this%we_lev_lon       )
    if (allocated(this%we_lev_lat       )) deallocate(this%we_lev_lat       )
    if (allocated(this%ptf_lon          )) deallocate(this%ptf_lon          )
    if (allocated(this%ptf_lat          )) deallocate(this%ptf_lat          )
    if (allocated(this%ptf_lev          )) deallocate(this%ptf_lev          )
    if (allocated(this%mfx_lon          )) deallocate(this%mfx_lon          )
    if (allocated(this%mfy_lat          )) deallocate(this%mfy_lat          )
    if (allocated(this%mfx_lat          )) deallocate(this%mfx_lat          )
    if (allocated(this%mfy_lon          )) deallocate(this%mfy_lon          )
    if (allocated(this%vor              )) deallocate(this%vor              )
    if (allocated(this%pv               )) deallocate(this%pv               )
    if (allocated(this%div              )) deallocate(this%div              )
    if (allocated(this%div2             )) deallocate(this%div2             )
    if (allocated(this%dmf              )) deallocate(this%dmf              )
    if (allocated(this%omg              )) deallocate(this%omg              )
    if (allocated(this%dudt_damp        )) deallocate(this%dudt_damp        )
    if (allocated(this%dvdt_damp        )) deallocate(this%dvdt_damp        )
    if (allocated(this%dptdt_damp       )) deallocate(this%dptdt_damp       )
    if (allocated(this%dqdt_damp        )) deallocate(this%dqdt_damp        )
    if (allocated(this%dudt_phys        )) deallocate(this%dudt_phys        )
    if (allocated(this%dvdt_phys        )) deallocate(this%dvdt_phys        )
    if (allocated(this%dptdt_phys       )) deallocate(this%dptdt_phys       )
    if (allocated(this%dqdt_phys        )) deallocate(this%dqdt_phys        )
    if (allocated(this%p_ptb            )) deallocate(this%p_ptb            )
    if (allocated(this%gz_ptb           )) deallocate(this%gz_ptb           )
    if (allocated(this%dp_ptb           )) deallocate(this%dp_ptb           )
    if (allocated(this%ad_ptb           )) deallocate(this%ad_ptb           )

  end subroutine aux_array_clear

  subroutine aux_array_final(this)

    type(aux_array_type), intent(inout) :: this

    call this%clear()

  end subroutine aux_array_final

end module dynamics_types_mod
