! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module dynamics_types_mod

  use const_mod
  use namelist_mod
  use latlon_field_types_mod
  use tracer_types_mod

  implicit none

  private

  public dstate_type
  public dtend_type
  public static_type
  public aux_array_type

  ! NOTE:
  !   Variables with '_lon', '_lat' and '_lev' are on the half grids on the corresponding direction.
  type dstate_type
    ! For nesting
    integer :: id = 0
    type(dstate_type), pointer :: parent => null()
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) u_lon
    type(latlon_field3d_type) v_lat
    type(latlon_field3d_type) we_lev
    type(latlon_field3d_type) gz
    type(latlon_field3d_type) gz_lev
    type(latlon_field3d_type) dmg
    type(latlon_field3d_type) dmg_lev
    type(latlon_field3d_type) pt
    type(latlon_field3d_type) t
    type(latlon_field3d_type) tv
    type(latlon_field3d_type) mg
    type(latlon_field3d_type) mg_lev
    type(latlon_field2d_type) mgs
    type(latlon_field3d_type) ph
    type(latlon_field3d_type) ph_lev
    type(latlon_field2d_type) phs
    type(latlon_field3d_type) rhod
    ! Nonhydrostatic variable
    type(latlon_field3d_type) we
    type(latlon_field3d_type) w
    type(latlon_field3d_type) w_lev
    type(latlon_field3d_type) p
    type(latlon_field3d_type) p_lev
    type(latlon_field2d_type) ps
    ! Total diagnostics
    real(r8) tm
    real(r8) te, te_ke, te_ie, te_pe
    real(r8) tpe
  contains
    procedure :: init  => dstate_init
    procedure :: clear => dstate_clear
    procedure :: c2a   => dstate_c2a
    procedure :: a2c   => dstate_a2c
    final dstate_final
  end type dstate_type

  type dtend_type
    type(latlon_field3d_type) du
    type(latlon_field3d_type) dv
    type(latlon_field3d_type) dgz
    type(latlon_field3d_type) dpt
    type(latlon_field2d_type) dmgs
#ifdef OUTPUT_H1_DTEND
    type(latlon_field3d_type) dudt_coriolis
    type(latlon_field3d_type) dvdt_coriolis
    type(latlon_field3d_type) dudt_wedudeta
    type(latlon_field3d_type) dvdt_wedvdeta
    type(latlon_field3d_type) dudt_dkedx
    type(latlon_field3d_type) dvdt_dkedy
    type(latlon_field3d_type) dudt_pgf
    type(latlon_field3d_type) dvdt_pgf
#endif
    logical :: update_u   = .false.
    logical :: update_v   = .false.
    logical :: update_gz  = .false.
    logical :: update_pt  = .false.
    logical :: update_mgs = .false.
    logical :: copy_gz    = .false.
    logical :: copy_pt    = .false.
    logical :: copy_mgs   = .false.
  contains
    procedure :: init        => dtend_init
    procedure :: reset_flags => dtend_reset_flags
    procedure :: clear       => dtend_clear
    final dtend_final
  end type dtend_type

  type static_type
    type(latlon_field2d_type) landmask
    ! Topography
    type(latlon_field2d_type) gzs
    type(latlon_field2d_type) zs_std
    type(latlon_field2d_type) dzsdx
    type(latlon_field2d_type) dzsdy
    ! Reference surface pressure
    type(latlon_field2d_type) ref_ps
    type(latlon_field2d_type) ref_ps_smth
    type(latlon_field2d_type) ref_ps_perb
    ! Coriolis parameters
    real(8), allocatable, dimension(:  ) :: f_lon
    real(8), allocatable, dimension(:  ) :: f_lat
    ! Weight for constructing tangential wind
    real(8), allocatable, dimension(:,:) :: tg_wgt_lon
    real(8), allocatable, dimension(:,:) :: tg_wgt_lat
  contains
    procedure :: init_stage1 => static_init_stage1
    procedure :: init_stage2 => static_init_stage2
    procedure :: clear       => static_clear
    final static_final
  end type static_type

  type aux_array_type
    ! Smagorinsky damping variables
    type(latlon_field3d_type) smag_t            ! tension strain
    type(latlon_field3d_type) smag_s            ! shear strain on vertex
    type(latlon_field3d_type) kmh               ! nonlinear diffusion coef
    type(latlon_field3d_type) kmh_lon           ! nonlinear diffusion coef on zonal edge
    type(latlon_field3d_type) kmh_lat           ! nonlinear diffusion coef on meridional edge
    ! Other variables
    type(latlon_field3d_type) v_lon             ! Meridional wind speed at lon edge (m s-1)
    type(latlon_field3d_type) u_lat             ! Zonal wind speed at lat edge (m s-1)
    type(latlon_field3d_type) ke                ! Kinetic energy
    type(latlon_field3d_type) pv_lon            ! Potential vorticity on zonal edge
    type(latlon_field3d_type) pv_lat            ! Potential vorticity on merdional edge
    type(latlon_field3d_type) dmg_lon           ! Mass on zonal edge
    type(latlon_field3d_type) dmg_lat           ! Mass on merdional edge
    type(latlon_field3d_type) dmg_vtx           ! Mass on vertex
    type(latlon_field3d_type) pkh_lev           ! Exner pressure on half levels
    type(latlon_field3d_type) we_lev_lon        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    type(latlon_field3d_type) we_lev_lat        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    type(latlon_field3d_type) ptf_lon           ! Potential temperature on the zonal edge
    type(latlon_field3d_type) ptf_lat           ! Potential temperature on the merdional edge
    type(latlon_field3d_type) ptf_lev           ! Potential temperature on the vertical edge
    type(latlon_field3d_type) mfx_lon           ! Normal mass flux on zonal edge
    type(latlon_field3d_type) mfy_lat           ! Normal mass flux on merdional edge
    type(latlon_field3d_type) mfx_lat           ! Tangient mass flux on zonal edge
    type(latlon_field3d_type) mfy_lon           ! Tangient mass flux on merdional edge
    type(latlon_field3d_type) vor               ! Vorticity (s-1)
    type(latlon_field3d_type) pv                ! Potential vorticity
    type(latlon_field3d_type) div               ! Divergence (s-1)
    type(latlon_field3d_type) div2              ! Laplacian of divergence (s-1)
    type(latlon_field3d_type) dmf               ! Mass flux divergence on full level (Pa s-1)
    type(latlon_field3d_type) omg               ! Vertical pressure velocity (Pa s-1)
    ! Tendencies from physics
    type(latlon_field3d_type) dudt_phys
    type(latlon_field3d_type) dvdt_phys
    type(latlon_field3d_type) dptdt_phys
    type(latlon_field4d_type) dqdt_phys
    ! Perturbed quantities for calculating HPGF
    type(latlon_field3d_type) p_ptb
    type(latlon_field3d_type) gz_ptb
    type(latlon_field3d_type) dp_ptb
    type(latlon_field3d_type) ad_ptb
    ! Nonhydrostatic variables
    type(latlon_field3d_type) u_lev_lon
    type(latlon_field3d_type) v_lev_lat
    type(latlon_field3d_type) mfx_lev_lon
    type(latlon_field3d_type) mfy_lev_lat
    type(latlon_field3d_type) adv_w_lev
    type(latlon_field3d_type) adv_gz_lev
  contains
    procedure :: init      => aux_array_init
    procedure :: init_phys => aux_array_init_phys
    procedure :: clear     => aux_array_clear
    final aux_array_final
  end type aux_array_type

contains

  subroutine dstate_init(this, filter_mesh, filter_halo, mesh, halo)

    class(dstate_type), intent(inout), target :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    name      = 'u'
    long_name = 'U wind component'
    units     = 'm s-1'
    call this%u%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'v'
    long_name = 'V wind component'
    units     = 'm s-1'
    call this%v%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'u_lon'
    long_name = 'U wind component on lon edge'
    units     = 'm s-1'
    call this%u_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'v_lat'
    long_name = 'V wind component on lat edge'
    units     = 'm s-1'
    call this%v_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'we_lev'
    long_name = 'Vertical coordinate velocity multiplied by dmg/deta on half level'
    units     = 'Pa s-1'
    if (baroclinic .or. advection) then
      call this%we_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'gz'
    long_name = 'Geopotential'
    units     = 'm2 s-2'
    call this%gz%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'gz_lev'
    long_name = 'Geopotential on half level'
    units     = 'm2 s-2'
    if (nonhydrostatic) then
      call this%gz_lev%init(name, long_name, units, 'lev', filter_mesh, filter_halo, halo_cross_pole=.true.)
    else
      call this%gz_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'dmg'
    long_name = 'Dry-air weight between two half levels'
    units     = 'Pa'
    call this%dmg%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'dmg_lev'
    long_name = 'Dry-air weight between two full levels'
    units     = 'Pa'
    if (baroclinic .or. advection) then
      call this%dmg_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 't'
    long_name = 'Temperature'
    units     = 'K'
    if (baroclinic) then
      call this%t%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'tv'
    long_name = 'Virtual temperature'
    units     = 'K'
    if (baroclinic .or. advection) then
      call this%tv%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'mg'
    long_name = 'Dry-air weight'
    units     = 'Pa'
    call this%mg%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'mg_lev'
    long_name = 'Dry-air weight on half level'
    units     = 'Pa'
    if (baroclinic .or. advection) then
      call this%mg_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'mgs'
    long_name = 'Dry-air weight on surface'
    units     = 'Pa'
    if (baroclinic .or. advection) then
      call this%mgs%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'ph'
    long_name = 'Hydrostatic pressure'
    units     = 'Pa'
    if (baroclinic .or. advection) then
      call this%ph%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'ph_lev'
    long_name = 'Hydrostatic pressure on half level'
    units     = 'Pa'
    if (baroclinic .or. advection) then
      call this%ph_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'rhod'
    long_name = 'Dry-air density'
    units     = 'kg m-3'
    if (baroclinic) then
      call this%rhod%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'pt'
    long_name = 'Modified potential temperature'
    units     = 'K'
    if (baroclinic) then
      call this%pt%init(name, long_name, units, 'cell', filter_mesh, filter_halo, halo_cross_pole=.true.)
    end if

    name      = 'phs'
    long_name = 'Hydrostatic pressure on surface'
    units     = 'Pa'
    if (baroclinic) then
      call this%phs%init(name, long_name, units, 'cell', mesh, halo)
      call this%phs%link(this%ph_lev, mesh%half_nlev)
    end if

    name      = 'we'
    long_name = 'Vertical coordinate velocity multiplied by dmg/deta on full level'
    units     = 'Pa s-1'
    if (nonhydrostatic) then
      call this%we%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'w'
    long_name = 'Vertical wind speed'
    units     = 'm s-1'
    if (nonhydrostatic) then
      call this%w%init(name, long_name, units, 'cell', mesh, halo)
    end if

    name      = 'w_lev'
    long_name = 'Vertical wind speed on half level'
    units     = 'm s-1'
    if (nonhydrostatic) then
      call this%w_lev%init(name, long_name, units, 'lev', filter_mesh, filter_halo, halo_cross_pole=.true.)
    end if

    name      = 'p'
    long_name = 'Pressure'
    units     = 'Pa'
    if (nonhydrostatic) then
      call this%p%init(name, long_name, units, 'cell', mesh, halo)
    else if (baroclinic) then
      call this%p%init(name, long_name, units, 'cell', mesh, halo, ptr_to=this%ph)
    end if

    name      = 'p_lev'
    long_name = 'Pressure on half level'
    units     = 'Pa'
    if (nonhydrostatic) then
      call this%p_lev%init(name, long_name, units, 'lev', mesh, halo)
    else if (baroclinic) then
      call this%p_lev%init(name, long_name, units, 'lev', mesh, halo, ptr_to=this%ph_lev)
    end if

    name      = 'ps'
    long_name = 'Surface pressure'
    units     = 'Pa'
    if (baroclinic) then
      call this%ps%init(name, long_name, units, 'cell', mesh, halo)
    end if

  end subroutine dstate_init

  subroutine dstate_clear(this)

    class(dstate_type), intent(inout) :: this

    call this%u      %clear()
    call this%v      %clear()
    call this%u_lon  %clear()
    call this%v_lat  %clear()
    call this%we_lev %clear()
    call this%gz     %clear()
    call this%gz_lev %clear()
    call this%dmg    %clear()
    call this%dmg_lev%clear()
    call this%pt     %clear()
    call this%t      %clear()
    call this%tv     %clear()
    call this%mg     %clear()
    call this%mg_lev %clear()
    call this%mgs    %clear()
    call this%ph     %clear()
    call this%ph_lev %clear()
    call this%phs    %clear()
    call this%rhod   %clear()
    call this%we     %clear()
    call this%w      %clear()
    call this%w_lev  %clear()
    call this%p      %clear()
    call this%p_lev  %clear()
    call this%ps     %clear()

  end subroutine dstate_clear

  subroutine dstate_final(this)

    type(dstate_type), intent(inout) :: this

    call this%clear()

  end subroutine dstate_final

  subroutine dstate_a2c(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    associate (mesh  => this%u%mesh, &
               u     => this%u     , & ! in
               v     => this%v     , & ! in
               u_lon => this%u_lon , & ! out
               v_lat => this%v_lat )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = 0.5_r8 * (u%d(i,j,k) + u%d(i+1,j,k))
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = 0.5_r8 * (v%d(i,j,k) + v%d(i,j+1,k))
        end do
      end do
    end do
    end associate

  end subroutine dstate_a2c

  subroutine dstate_c2a(this)

    class(dstate_type), intent(inout) :: this

    integer i, j, k

    associate (mesh  => this%u%mesh, &
               u     => this%u     , & ! out
               v     => this%v     , & ! out
               u_lon => this%u_lon , & ! in
               v_lat => this%v_lat )   ! in
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          u%d(i,j,k) = 0.5_r8 * (u_lon%d(i,j,k) + u_lon%d(i-1,j,k))
          v%d(i,j,k) = 0.5_r8 * (v_lat%d(i,j,k) + v_lat%d(i,j-1,k))
        end do
      end do
    end do
    end associate

  end subroutine dstate_c2a

  subroutine dtend_init(this, filter_mesh, filter_halo, mesh, halo)

    class(dtend_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    name      = 'dudt'
    long_name = 'Dynamic tendency of U wind component'
    units     = 'm s-2'
    call this%du%init(name, long_name, units, 'lon', filter_mesh, filter_halo)

    name      = 'dvdt'
    long_name = 'Dynamic tendency of V wind component'
    units     = 'm s-2'
    call this%dv%init(name, long_name, units, 'lat', filter_mesh, filter_halo)

    name      = 'dptdt'
    long_name = 'Dynamic tendency of modified potential temperature'
    units     = 'K s-1'
    if (baroclinic) then
      call this%dpt%init(name, long_name, units, 'cell', filter_mesh, filter_halo)
    end if

    name      = 'dmgsdt'
    long_name = 'Dynamic tendency of dry-air weight on surface'
    units     = 'Pa s-1'
    if (baroclinic) then
      call this%dmgs%init(name, long_name, units, 'cell', filter_mesh, filter_halo)
    end if

    name      = 'dgzdt'
    long_name = 'Dynamic tendency of geopotential'
    units     = 'm2 s-2'
    if (nonhydrostatic .or. .not. baroclinic) then
      call this%dgz%init(name, long_name, units, 'cell', filter_mesh, filter_halo)
    end if

#ifdef OUTPUT_H1_DTEND
    name      = 'dudt_coriolis'
    long_name = 'Dynamic tendency of U wind component due to Coriolis force'
    units     = 'm s-2'
    call this%dudt_coriolis%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'dvdt_coriolis'
    long_name = 'Dynamic tendency of V wind component due to Coriolis force'
    units     = 'm s-2'
    call this%dvdt_coriolis%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'dudt_wedudeta'
    long_name = 'Dynamic tendency of U wind component due to vertical advection'
    units     = 'm s-2'
    if (baroclinic) then
      call this%dudt_wedudeta%init(name, long_name, units, 'lon', mesh, halo)
    end if

    name      = 'dvdt_wedvdeta'
    long_name = 'Dynamic tendency of V wind component due to vertical advection'
    units     = 'm s-2'
    if (baroclinic) then
      call this%dvdt_wedvdeta%init(name, long_name, units, 'lat', mesh, halo)
    end if

    name      = 'dudt_dkedx'
    long_name = 'Dynamic tendency of U wind component due to kinetic gradient'
    units     = 'm s-2'
    call this%dudt_dkedx%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'dvdt_dkedy'
    long_name = 'Dynamic tendency of V wind component due to kinetic gradient'
    units     = 'm s-2'
    call this%dvdt_dkedy%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'dudt_pgf'
    long_name = 'Dynamic tendency of U wind component due to pressure gradient force'
    units     = 'm s-2'
    call this%dudt_pgf%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'dvdt_pgf'
    long_name = 'Dynamic tendency of V wind component due to pressure gradient force'
    units     = 'm s-2'
    call this%dvdt_pgf%init(name, long_name, units, 'lat', mesh, halo)
#endif

  end subroutine dtend_init

  subroutine dtend_reset_flags(this)

    class(dtend_type), intent(inout) :: this

    this%du%d = 0
    this%dv%d = 0

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

    call this%du  %clear()
    call this%dv  %clear()
    call this%dgz %clear()
    call this%dpt %clear()
    call this%dmgs%clear()

#ifdef OUTPUT_H1_DTEND
    call this%dudt_coriolis%clear()
    call this%dvdt_coriolis%clear()
    call this%dudt_wedudeta%clear()
    call this%dvdt_wedvdeta%clear()
    call this%dudt_dkedx   %clear()
    call this%dvdt_dkedy   %clear()
    call this%dudt_pgf     %clear()
    call this%dvdt_pgf     %clear()
#endif

  end subroutine dtend_clear

  subroutine dtend_final(this)

    type(dtend_type), intent(inout) :: this

    call this%clear()

  end subroutine dtend_final

  subroutine static_init_stage1(this, filter_mesh, filter_halo, mesh, halo)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    name      = 'gzs'
    long_name = 'Surface geopotential'
    units     = 'm2 s-2'
    call this%gzs%init(name, long_name, units, 'cell', filter_mesh, filter_halo)

    name      = 'landmask'
    long_name = 'Land mask'
    units     = '1'
    call this%landmask%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'zs_std'
    long_name = 'Subgrid variance of surface geopotential height'
    units     = 'm2 s-2'
    call this%zs_std%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'dzsdx'
    long_name = 'Zonal gradient of surface geopotential height'
    units     = '1'
    call this%dzsdx%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'dzsdy'
    long_name = 'Meridional gradient of surface geopotential height'
    units     = '1'
    call this%dzsdy%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'ref_ps'
    long_name = 'Reference surface pressure'
    units     = 'Pa'
    call this%ref_ps%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_ps_smth'
    long_name = 'Smoothed reference surface pressure'
    units     = 'Pa'
    call this%ref_ps_smth%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_ps_perb'
    long_name = 'Perturbation of reference surface pressure'
    units     = 'Pa'
    call this%ref_ps_perb%init(name, long_name, units, 'cell', mesh, halo)

    allocate(this%f_lon       (mesh%full_jms:mesh%full_jme)); this%f_lon      = inf
    allocate(this%f_lat       (mesh%half_jms:mesh%half_jme)); this%f_lat      = inf
    allocate(this%tg_wgt_lon(2,mesh%full_jms:mesh%full_jme)); this%tg_wgt_lon = inf
    allocate(this%tg_wgt_lat(2,mesh%half_jms:mesh%half_jme)); this%tg_wgt_lat = inf

  end subroutine static_init_stage1

  subroutine static_init_stage2(this, mesh)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: mesh

    integer j

    do j = mesh%full_jds, mesh%full_jde
      this%f_lon(j) = 2 * omega * mesh%full_sin_lat(j)
    end do
    do j = mesh%half_jds, mesh%half_jde
      this%f_lat(j) = 2 * omega * mesh%half_sin_lat(j)
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
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        this%tg_wgt_lon(1,j) = mesh%le_lat(j-1) / mesh%de_lon(j) * 0.25d0
        this%tg_wgt_lon(2,j) = mesh%le_lat(j  ) / mesh%de_lon(j) * 0.25d0
      end do

      do j = mesh%half_jds, mesh%half_jde
        this%tg_wgt_lat(1,j) = mesh%le_lon(j  ) / mesh%de_lat(j) * 0.25d0
        this%tg_wgt_lat(2,j) = mesh%le_lon(j+1) / mesh%de_lat(j) * 0.25d0
      end do
    case ('th09')
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        this%tg_wgt_lon(1,j) = mesh%le_lat(j-1) / mesh%de_lon(j) * mesh%area_subcell(2,j  ) / mesh%area_cell(j  )
        this%tg_wgt_lon(2,j) = mesh%le_lat(j  ) / mesh%de_lon(j) * mesh%area_subcell(1,j  ) / mesh%area_cell(j  )
      end do

      do j = mesh%half_jds, mesh%half_jde
        this%tg_wgt_lat(1,j) = mesh%le_lon(j  ) / mesh%de_lat(j) * mesh%area_subcell(1,j  ) / mesh%area_cell(j  )
        this%tg_wgt_lat(2,j) = mesh%le_lon(j+1) / mesh%de_lat(j) * mesh%area_subcell(2,j+1) / mesh%area_cell(j+1)
      end do
    end select

  end subroutine static_init_stage2

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    call this%gzs        %clear()
    call this%landmask   %clear()
    call this%zs_std     %clear()
    call this%dzsdx      %clear()
    call this%dzsdy      %clear()
    call this%ref_ps     %clear()
    call this%ref_ps_smth%clear()
    call this%ref_ps_perb%clear()

    if (allocated(this%f_lon     )) deallocate(this%f_lon     )
    if (allocated(this%f_lat     )) deallocate(this%f_lat     )
    if (allocated(this%tg_wgt_lon)) deallocate(this%tg_wgt_lon)
    if (allocated(this%tg_wgt_lat)) deallocate(this%tg_wgt_lat)

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

  subroutine aux_array_init(this, filter_mesh, filter_halo, mesh, halo)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    if (use_smag_damp) then
      name      = 'smag_t'
      long_name = 'Tension of horizontal wind for Smagorinsky damping'
      units     = 's-2'
      call this%smag_t%init(name, long_name, units, 'cell', mesh, halo)

      name      = 'smag_s'
      long_name = 'Shear of horizontal wind for Smagorinsky damping'
      units     = 's-2'
      call this%smag_s%init(name, long_name, units, 'vtx', mesh, halo)

      name      = 'kmh'
      long_name = 'Horizontal eddy viscosity for Smagorinsky damping'
      units     = 's-1'
      call this%kmh%init(name, long_name, units, 'cell', mesh, halo)

      name      = 'kmh_lon'
      long_name = 'Horizontal eddy viscosity for Smagorinsky damping on lon edge'
      units     = 's-1'
      call this%kmh_lon%init(name, long_name, units, 'lon', mesh, halo)

      name      = 'kmh_lat'
      long_name = 'Horizontal eddy viscosity for Smagorinsky damping on lat edge'
      units     = 's-1'
      call this%kmh_lat%init(name, long_name, units, 'lat', mesh, halo)
    end if

    name      = 'v_lon'
    long_name = 'V wind component on lon edge'
    units     = 'm s-1'
    call this%v_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'u_lat'
    long_name = 'U wind component on lat edge'
    units     = 'm s-1'
    call this%u_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'ke'
    long_name = 'Kinetic energy'
    units     = 'm2 s-2'
    call this%ke%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'pv_lon'
    long_name = 'Potential vorticity on lon edge'
    units     = 'Pa-1 s-1'
    call this%pv_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'pv_lat'
    long_name = 'Potential vorticity on lat edge'
    units     = 'Pa-1 s-1'
    call this%pv_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'dmg_lon'
    long_name = 'Dry-air weight between two half levels on lon edge'
    units     = 'Pa'
    call this%dmg_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'dmg_lat'
    long_name = 'Dry-air weight between two half levels on lat edge'
    units     = 'Pa'
    call this%dmg_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'dmg_vtx'
    long_name = 'Dry-air weight between two half levels on vtx edge'
    units     = 'Pa'
    call this%dmg_vtx%init(name, long_name, units, 'vtx', mesh, halo)

    name      = 'pkh_lev'
    long_name = 'Hydrostatic pressure under Kappa exponent on half level'
    units     = 'Pa'
    if (baroclinic) then
      call this%pkh_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'we_lev_lon'
    long_name = 'Vertical coordinate velocity multiplied by dmg/deta on lon edge'
    units     = 'Pa s-1'
    if (baroclinic) then
      call this%we_lev_lon%init(name, long_name, units, 'lev_lon', mesh, halo)
    end if

    name      = 'we_lev_lat'
    long_name = 'Vertical coordinate velocity multiplied by dmg/deta on lat edge'
    units     = 'Pa s-1'
    if (baroclinic) then
      call this%we_lev_lat%init(name, long_name, units, 'lev_lat', mesh, halo)
    end if

    name      = 'ptf_lon'
    long_name = 'Modified potential temperature flux on lon edge'
    units     = 'K m s-1'
    if (baroclinic) then
      call this%ptf_lon%init(name, long_name, units, 'lon', mesh, halo)
    end if

    name      = 'ptf_lat'
    long_name = 'Modified potential temperature flux on lat edge'
    units     = 'K m s-1'
    if (baroclinic) then
      call this%ptf_lat%init(name, long_name, units, 'lat', mesh, halo)
    end if

    name      = 'ptf_lev'
    long_name = 'Modified potential temperature flux on half level'
    units     = 'K m s-1'
    if (baroclinic) then
      call this%ptf_lev%init(name, long_name, units, 'lev', mesh, halo)
    end if

    name      = 'mfx_lon'
    long_name = 'Zonal mass flux on lon edge'
    units     = 'Pa m s-1'
    call this%mfx_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'mfy_lat'
    long_name = 'Meridional mass flux on lat edge'
    units     = 'Pa m s-1'
    call this%mfy_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'mfx_lat'
    long_name = 'Zonal mass flux on lat edge'
    units     = 'Pa m s-1'
    call this%mfx_lat%init(name, long_name, units, 'lat', mesh, halo)

    name      = 'mfy_lon'
    long_name = 'Meridional mass flux on lon edge'
    units     = 'Pa m s-1'
    call this%mfy_lon%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'vor'
    long_name = 'Relative vorticity'
    units     = 's-1'
    call this%vor%init(name, long_name, units, 'vtx', mesh, halo)

    name      = 'pv'
    long_name = 'Potential vorticity'
    units     = 'Pa-1 s-1'
    call this%pv%init(name, long_name, units, 'vtx', mesh, halo, halo_cross_pole=.true.)

    name      = 'div'
    long_name = 'Divergence'
    units     = 's-1'
    call this%div%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'div2'
    long_name = 'Gradient of divergence'
    units     = 'm-1 s-1'
    call this%div2%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'dmf'
    long_name = 'Mass flux divergence'
    units     = 'Pa s-1'
    call this%dmf%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'omg'
    long_name = 'Omega'
    units     = 'Pa s-1'
    if (baroclinic) then
      call this%omg%init(name, long_name, units, 'cell', mesh, halo)
    end if

    if (pgf_scheme == 'ptb') then
      name      = 'p_ptb'
      long_name = 'Perturbation of pressure'
      units     = 'Pa'
      call this%p_ptb%init(name, long_name, units, 'cell', mesh, halo)

      name      = 'gz_ptb'
      long_name = 'Perturbation of geopotential'
      units     = 'm2 s-2'
      call this%gz_ptb%init(name, long_name, units, 'cell', mesh, halo)

      name      = 'dp_ptb'
      long_name = 'Perturbation of dry-air weight'
      units     = 'Pa'
      call this%dp_ptb%init(name, long_name, units, 'cell', mesh, halo)

      name      = 'ad_ptb'
      long_name = 'Perturbation of specific density of dry-air'
      units     = 'kg-1 m3'
      call this%ad_ptb%init(name, long_name, units, 'cell', mesh, halo)
    end if

    if (nonhydrostatic) then
      name      = 'u_lev_lon'
      long_name = 'U wind component on lon edge on half level'
      units     = 'm s-1'
      call this%u_lev_lon%init(name, long_name, units, 'lev_lon', mesh, halo)

      name      = 'v_lev_lat'
      long_name = 'V wind component on lat edge on half level'
      units     = 'm s-1'
      call this%v_lev_lat%init(name, long_name, units, 'lev_lat', mesh, halo)

      name      = 'mfx_lev_lon'
      long_name = 'Zonal mass flux on lon edge on half level'
      units     = 'Pa m s-1'
      call this%mfx_lev_lon%init(name, long_name, units, 'lev_lon', mesh, halo)

      name      = 'mfy_lev_lat'
      long_name = 'Meridional mass flux on lat edge on half level'
      units     = 'Pa m s-1'
      call this%mfy_lev_lat%init(name, long_name, units, 'lev_lat', mesh, halo)

      name      = 'adv_w_lev'
      long_name = 'Advection tendency of vertical wind speed on half level'
      units     = 'm s-2'
      call this%adv_w_lev%init(name, long_name, units, 'lev', filter_mesh, filter_halo)

      name      = 'adv_gz_lev'
      long_name = 'Advection tendency of geopotential on half level'
      units     = 'm2 s-2'
      call this%adv_gz_lev%init(name, long_name, units, 'lev', filter_mesh, filter_halo)
    end if

  end subroutine aux_array_init

  subroutine aux_array_init_phys(this, filter_mesh, filter_halo, mesh, halo)

    class(aux_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    if (trim(physics_suite) /= '') then
      name      = 'dudt_phys'
      long_name = 'Physics tendency of U wind component'
      units     = 'm s-2'
      if (filter_ptend) then
        call this%dudt_phys%init(name, long_name, units, 'lon', filter_mesh, filter_halo)
      else
        call this%dudt_phys%init(name, long_name, units, 'lon', mesh, halo)
      end if

      name      = 'dvdt_phys'
      long_name = 'Physics tendency of V wind component'
      units     = 'm s-2'
      if (filter_ptend) then
        call this%dvdt_phys%init(name, long_name, units, 'lat', filter_mesh, filter_halo)
      else
        call this%dvdt_phys%init(name, long_name, units, 'lat', mesh, halo)
      end if

      name      = 'dptdt_phys'
      long_name = 'Physics tendency of modified potential temperature'
      units     = 'K s-1'
      if (filter_ptend) then
        call this%dptdt_phys%init(name, long_name, units, 'cell', filter_mesh, filter_halo)
      else
        call this%dptdt_phys%init(name, long_name, units, 'cell', mesh, halo)
      end if

      name      = 'dqdt_phys'
      long_name = 'Physics tendency of tracer dry mixing ratio'
      units     = 'kg kg-1 s-1'
      if (filter_ptend) then
        call this%dqdt_phys%init(name, long_name, units, 'cell', filter_mesh, filter_halo, n4=ntracers)
      else
        call this%dqdt_phys%init(name, long_name, units, 'cell', mesh, halo, n4=ntracers)
      end if
    end if

  end subroutine aux_array_init_phys

  subroutine aux_array_clear(this)

    class(aux_array_type), intent(inout) :: this

    call this%smag_t     %clear()
    call this%smag_s     %clear()
    call this%kmh        %clear()
    call this%kmh_lon    %clear()
    call this%kmh_lat    %clear()
    call this%v_lon      %clear()
    call this%u_lat      %clear()
    call this%ke         %clear()
    call this%pv_lon     %clear()
    call this%pv_lat     %clear()
    call this%dmg_lon    %clear()
    call this%dmg_lat    %clear()
    call this%dmg_vtx    %clear()
    call this%pkh_lev    %clear()
    call this%we_lev_lon %clear()
    call this%we_lev_lat %clear()
    call this%ptf_lon    %clear()
    call this%ptf_lat    %clear()
    call this%ptf_lev    %clear()
    call this%mfx_lon    %clear()
    call this%mfy_lat    %clear()
    call this%mfx_lat    %clear()
    call this%mfy_lon    %clear()
    call this%vor        %clear()
    call this%pv         %clear()
    call this%div        %clear()
    call this%div2       %clear()
    call this%dmf        %clear()
    call this%omg        %clear()
    call this%dudt_phys  %clear()
    call this%dvdt_phys  %clear()
    call this%dptdt_phys %clear()
    call this%dqdt_phys  %clear()
    call this%p_ptb      %clear()
    call this%gz_ptb     %clear()
    call this%dp_ptb     %clear()
    call this%ad_ptb     %clear()
    call this%u_lev_lon  %clear()
    call this%v_lev_lat  %clear()
    call this%mfx_lev_lon%clear()
    call this%mfy_lev_lat%clear()
    call this%adv_w_lev  %clear()
    call this%adv_gz_lev %clear()

  end subroutine aux_array_clear

  subroutine aux_array_final(this)

    type(aux_array_type), intent(inout) :: this

    call this%clear()

  end subroutine aux_array_final

end module dynamics_types_mod
