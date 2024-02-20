! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module namelist_mod

  use string
  use flogger
  use const_mod, only: r8, const_init, time_scale, earth_day_seconds, mars_sol_seconds

  implicit none

  character(30)   :: planet               = 'earth'
  logical         :: use_aqua_planet      = .false.

  integer         :: start_time(5)        = 0
  integer         :: end_time(5)          = 0
  real(r8)        :: run_hours            = 0
  real(r8)        :: run_days             = 0
  real(r8)        :: run_years            = 0
  real(r8)        :: run_my               = 0
  real(r8)        :: run_sol              = 0

  ! Individual time step sizes in target planet time system
  real(r8)        :: dt_dyn               = 0
  real(r8)        :: dt_adv               = 0
  real(r8)        :: dt_phys              = 0
  ! Physics-dynamics coupling type
  ! 1: Update dynamics and tracers after their own calculation.
  ! 2: Update dynamics and tracers after advection (moisture).
  ! 3: Update dynamics and tracers after physics.
  ! 4: Update dynamics at RK sub-steps and tracers after advection.
  integer         :: pdc_type             = 2

  character(256)  :: case_desc            = 'N/A'
  character(256)  :: case_name            = 'N/A'
  character(30 )  :: test_case            = 'N/A'
  character(30)   :: initial_time         = '2000-01-01T00:00:00'
  character(30 )  :: history_interval(1)  = 'N/A'
  character(30 )  :: restart_interval     = 'N/A'
  character(30 )  :: print_interval       = '1 hours'
  character(256)  :: initial_file         = 'N/A'
  character(256)  :: restart_file         = 'N/A'
  character(256)  :: topo_file            = 'N/A'
  character(30 )  :: topo_type            = 'etopo1' ! etopo1, gmted, mola32
  character(256)  :: bkg_file             = 'N/A'
  character(30 )  :: bkg_type             = 'era5'

  integer         :: nlon
  integer         :: nlat
  integer         :: nlev                 = 1

  logical         :: baroclinic           = .false.
  logical         :: hydrostatic          = .true.
  logical         :: nonhydrostatic       = .false.
  logical         :: advection            = .false.
  logical         :: restart              = .false.

  character(30)   :: physics_suite        = 'N/A'
  character(30)   :: mp_scheme            = 'N/A'
  character(30)   :: pbl_scheme           = 'N/A'
  character(256)  :: cam_namelist_path    = 'N/A'
  logical         :: filter_ptend         = .false.

  character(256)  :: gmcore_data_dir      = 'N/A'

  integer         :: nproc_x(20)          = 0
  integer         :: nproc_y(20)          = 0
  integer         :: lon_hw               = 3
  integer         :: lat_hw               = 3
  character(30)   :: proc_layout          = 'lon>lat' ! or 'lat>lon'

  character(30)   :: tangent_wgt_scheme   = 'classic'

  real(r8)        :: implicit_w_wgt       = 0.55_r8

  character(30)   :: vert_coord_scheme    = 'hybrid'
  character(30)   :: vert_coord_template  = 'N/A'
  character(30)   :: refer_state_scheme   = 'wrf'
  real(r8)        :: ptop                 = 2.194e2_r8
  real(r8)        :: hybrid_coord_p0      = 1.0e5_r8

  ! Parameters for generating hybrid levels from WRF.
  real(r8)        :: tiso                 = 300.0_r8  ! Isothermal temperature (K)
  real(r8)        :: dzbot                = 10.0_r8   ! Bottom layer thickness (m)
  real(r8)        :: dzmax                = 1000.0_r8 ! Maximum layer thickness (m)
  real(r8)        :: dzstretch_s          = 1.3_r8    ! Stretching factor for surface layers
  real(r8)        :: dzstretch_u          = 1.1_r8    ! Stretching factor for upper layers

  ! Parameters for generating hybrid levels from NCEP.
  real(r8)        :: hybrid_coord_ncep_psig   = 0
  real(r8)        :: hybrid_coord_ncep_ppre   = 0
  real(r8)        :: hybrid_coord_ncep_dpbot  = 0
  real(r8)        :: hybrid_coord_ncep_dpsig  = 0
  real(r8)        :: hybrid_coord_ncep_dppre  = 0
  real(r8)        :: hybrid_coord_ncep_dptop  = 0

  integer         :: ke_scheme            = 2
  real(r8)        :: ke_cell_wgt          = 0.5_r8

  character(30)   :: pv_scheme            = 'upwind' ! midpoint, upwind, ffsl
  logical         :: pv_pole_stokes       = .true.
  integer         :: upwind_order_pv      = 3
  real(r8)        :: upwind_wgt_pv        = 1

  character(8)    :: pgf_scheme           = 'lin97'
  integer         :: coriolis_scheme      = 1

  character(8)    :: pt_adv_scheme        = 'ffsl'
  character(8)    :: nh_adv_scheme        = 'upwind'
  character(8)    :: limiter_type         = 'mono'
  character(8)    :: ffsl_flux_type       = 'ppm'
  character(8)    :: tvd_limiter_type     = 'van_leer'

  character(8)    :: zonal_tridiag_solver = 'spk' ! mkl, spk

  integer         :: weno_order           = -1 ! -1, 3
  integer         :: upwind_order         = 3  ! -1, 1, 3
  real(r8)        :: upwind_wgt           = 1.0_r8
  real(r8)        :: upwind_wgt_pt        = 0.25_r8

  integer         :: vert_weno_order      = -1 ! -1, 3
  integer         :: vert_upwind_order    = 3  ! -1, 1, 3
  real(r8)        :: vert_upwind_wgt      = 1.0_r8

  character(30)   :: time_scheme          = 'wrfrk3'

  real(r8)        :: coarse_pole_mul      = 0
  real(r8)        :: coarse_pole_decay    = 100.0

  ! Filter settings
  real(r8)        :: max_wave_speed       = 300
  real(r8)        :: max_cfl              = 0.5
  real(r8)        :: filter_coef_a        = 1.5
  real(r8)        :: filter_coef_b        = 0.2
  real(r8)        :: filter_coef_c        = 0.5
  real(r8)        :: filter_min_width     = 0.0

  ! Damping settings
  logical         :: use_topo_smooth      = .false.
  real(r8)        :: topo_max_slope       = 0.12_r8
  integer         :: topo_smooth_cycles   = 3
  logical         :: use_div_damp         = .false.
  integer         :: div_damp_cycles      = 1
  integer         :: div_damp_order       = 2
  real(r8)        :: div_damp_top         = 3
  integer         :: div_damp_k0          = 6
  real(r8)        :: div_damp_pole        = 0
  real(r8)        :: div_damp_pole_x      = 10
  real(r8)        :: div_damp_pole_y      = 10
  real(r8)        :: div_damp_lat0        = 70
  real(r8)        :: div_damp_coef2       = 1.0_r8 / 128.0_r8
  real(r8)        :: div_damp_coef4       = 0.001_r8
  logical         :: use_vor_damp         = .false.
  integer         :: vor_damp_cycles      = 1
  integer         :: vor_damp_order       = 2
  real(r8)        :: vor_damp_coef2       = 0.001_r8
  real(r8)        :: vor_damp_top         = 1
  integer         :: vor_damp_k0          = 6
  real(r8)        :: vor_damp_pole        = 0
  real(r8)        :: vor_damp_pole_x      = 200
  real(r8)        :: vor_damp_pole_y      = 200
  real(r8)        :: vor_damp_lat0        = 60
  real(r8)        :: rayleigh_damp_w_coef = 0.2
  real(r8)        :: rayleigh_damp_top    = 10.0d3 ! m
  logical         :: use_smag_damp        = .false.
  integer         :: smag_damp_cycles     = 1
  real(r8)        :: smag_damp_coef       = 0.015

  ! Input settings
  integer         :: input_ngroup         = 0

  ! Output settings
#if (REAL_KIND == 4)
  character(8)    :: output_i0_dtype      = 'r4'
#elif (REAL_KIND == 8)
  character(8)    :: output_i0_dtype      = 'r8'
#endif
  logical         :: output_h0            = .true.
  character(8)    :: output_h0_dtype      = 'r4'
  logical         :: output_h1            = .false.
  character(30)   :: output_h0_new_file   = ''
  character(8)    :: output_h0_vars(100)  = ''
  integer         :: output_ngroup        = 0

  namelist /gmcore_control/     &
    planet                    , &
    use_aqua_planet           , &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    nlon                      , &
    nlat                      , &
    nlev                      , &
    nonhydrostatic            , &
    advection                 , &
    nproc_x                   , &
    nproc_y                   , &
    lon_hw                    , &
    lat_hw                    , &
    proc_layout               , &
    initial_time              , &
    start_time                , &
    end_time                  , &
    dt_dyn                    , &
    dt_adv                    , &
    dt_phys                   , &
    pdc_type                  , &
    run_years                 , &
    run_my                    , &
    run_sol                   , &
    run_hours                 , &
    run_days                  , &
    history_interval          , &
    restart_interval          , &
    print_interval            , &
    initial_file              , &
    output_i0_dtype           , &
    restart_file              , &
    restart                   , &
    topo_file                 , &
    topo_type                 , &
    bkg_file                  , &
    bkg_type                  , &
    tangent_wgt_scheme        , &
    implicit_w_wgt            , &
    vert_coord_scheme         , &
    vert_coord_template       , &
    refer_state_scheme        , &
    ptop                      , &
    hybrid_coord_p0           , &
    tiso                      , &
    dzbot                     , &
    dzmax                     , &
    dzstretch_s               , &
    dzstretch_u               , &
    hybrid_coord_ncep_psig    , &
    hybrid_coord_ncep_ppre    , &
    hybrid_coord_ncep_dpbot   , &
    hybrid_coord_ncep_dpsig   , &
    hybrid_coord_ncep_dppre   , &
    hybrid_coord_ncep_dptop   , &
    ke_scheme                 , &
    ke_cell_wgt               , &
    pv_scheme                 , &
    pv_pole_stokes            , &
    upwind_order_pv           , &
    upwind_wgt_pv             , &
    pgf_scheme                , &
    coriolis_scheme           , &
    pt_adv_scheme             , &
    nh_adv_scheme             , &
    limiter_type              , &
    ffsl_flux_type            , &
    tvd_limiter_type          , &
    zonal_tridiag_solver      , &
    weno_order                , &
    upwind_order              , &
    upwind_wgt                , &
    upwind_wgt_pt             , &
    vert_weno_order           , &
    vert_upwind_order         , &
    vert_upwind_wgt           , &
    time_scheme               , &
    max_wave_speed            , &
    max_cfl                   , &
    filter_coef_a             , &
    filter_coef_b             , &
    filter_coef_c             , &
    filter_min_width          , &
    coarse_pole_mul           , &
    coarse_pole_decay         , &
    physics_suite             , &
    mp_scheme                 , &
    pbl_scheme                , &
    cam_namelist_path         , &
    filter_ptend              , &
    gmcore_data_dir           , &
    use_topo_smooth           , &
    topo_max_slope            , &
    topo_smooth_cycles        , &
    use_div_damp              , &
    div_damp_cycles           , &
    div_damp_order            , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    div_damp_k0               , &
    div_damp_top              , &
    div_damp_pole             , &
    div_damp_pole_x           , &
    div_damp_pole_y           , &
    div_damp_lat0             , &
    use_vor_damp              , &
    vor_damp_cycles           , &
    vor_damp_order            , &
    vor_damp_k0               , &
    vor_damp_coef2            , &
    vor_damp_top              , &
    vor_damp_pole             , &
    vor_damp_pole_x           , &
    vor_damp_pole_y           , &
    vor_damp_lat0             , &
    rayleigh_damp_w_coef      , &
    rayleigh_damp_top         , &
    use_smag_damp             , &
    smag_damp_cycles          , &
    smag_damp_coef            , &
    input_ngroup              , &
    output_h0                 , &
    output_h0_dtype           , &
    output_h1                 , &
    output_h0_new_file        , &
    output_h0_vars            , &
    output_ngroup

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_control)
    close(10)

    ! Here we set baroclinic according to levels.
    baroclinic = nlev > 1
    if (.not. baroclinic) then
      hydrostatic    = .false.
      nonhydrostatic = .false.
      ke_scheme      = 1
      if (use_div_damp) then
        use_div_damp   = .false.
        call log_warning('In shallow water mode, no need to use divergence damping! Turn it off for you.')
      end if
    else
      hydrostatic = .not. nonhydrostatic
    end if

    if (advection) then
      hydrostatic    = .false.
      baroclinic     = .false.
      nonhydrostatic = .false.
    end if

    if (dt_dyn  == 0) dt_dyn  = dt_adv
    if (dt_adv  == 0) dt_adv  = dt_dyn
    if (dt_phys == 0) dt_phys = dt_adv

    ! Convert time step sizes to Earth time system
    select case (planet)
    case ('mars')
      time_scale = mars_sol_seconds / earth_day_seconds
    end select
    dt_dyn  = dt_dyn  * time_scale
    dt_adv  = dt_adv  * time_scale
    dt_phys = dt_phys * time_scale

    if (physics_suite == 'N/A') then
      pdc_type       = 0
    end if

    if (.not. use_div_damp) then
      div_damp_order = 0
    else if (div_damp_pole /= 0) then
      div_damp_pole_x = div_damp_pole
      div_damp_pole_y = div_damp_pole
    end if

    if (.not. use_vor_damp) then
      vor_damp_order = 0
    else if (vor_damp_pole /= 0) then
      vor_damp_pole_x = vor_damp_pole
      vor_damp_pole_y = vor_damp_pole
    end if

  end subroutine parse_namelist

  subroutine print_namelist()

      write(*, *) '=================== GMCORE Parameters ==================='
      write(*, *) 'case_name           = ', trim(case_name)
      write(*, *) 'nlon                = ', to_str(nlon)
      write(*, *) 'nlat                = ', to_str(nlat)
      write(*, *) 'nlev                = ', to_str(nlev)
    if (coarse_pole_mul /= 0) then
      write(*, *) 'coarse_pole_mul     = ', to_str(coarse_pole_mul, 3)
      write(*, *) 'coarse_pole_decay   = ', to_str(coarse_pole_decay, 3)
    end if
      write(*, *) 'physics_suite       = ', trim(physics_suite)
      write(*, *) 'mp_scheme           = ', trim(mp_scheme)
      write(*, *) 'pbl_scheme          = ', trim(pbl_scheme)
      write(*, *) 'hydrostatic         = ', to_str(hydrostatic)
      write(*, *) 'nonhydrostatic      = ', to_str(nonhydrostatic)
      write(*, *) 'vert_coord_scheme   = ', trim(vert_coord_scheme)
      write(*, *) 'vert_coord_template = ', trim(vert_coord_template)
      write(*, *) 'ptop                = ', to_str(ptop, 4)
      write(*, *) 'hybrid_coord_p0     = ', to_str(hybrid_coord_p0, 2)
      write(*, *) 'dt_dyn              = ', to_str(dt_dyn , 2)
      write(*, *) 'dt_adv              = ', to_str(dt_adv , 2)
      write(*, *) 'dt_phys             = ', to_str(dt_phys, 2)
      write(*, *) 'pdc_type            = ', to_str(pdc_type)
      write(*, *) 'max_wave_speed      = ', to_str(max_wave_speed, 2)
      write(*, *) 'max_cfl             = ', to_str(max_cfl, 2)
      write(*, *) 'filter_coef_a       = ', filter_coef_a
      write(*, *) 'filter_coef_b       = ', filter_coef_b
      write(*, *) 'filter_coef_c       = ', filter_coef_c
      write(*, *) 'filter_min_width    = ', filter_min_width
      write(*, *) 'filter_ptend        = ', to_str(filter_ptend)
      write(*, *) 'pgf_scheme          = ', trim(pgf_scheme)
      write(*, *) 'pt_adv_scheme       = ', trim(pt_adv_scheme)
      write(*, *) 'nh_adv_scheme       = ', trim(nh_adv_scheme)
      write(*, *) 'limiter_type        = ', trim(limiter_type)
    if (pt_adv_scheme == 'ffsl') then
      write(*, *) 'ffsl_flux_type      = ', trim(ffsl_flux_type)
    end if
      write(*, *) 'ke_scheme           = ', to_str(ke_scheme)
    if (ke_scheme == 2) then
      write(*, *) 'ke_cell_wgt         = ', to_str(ke_cell_wgt, 2)
    end if
      write(*, *) 'pv_scheme           = ', trim(pv_scheme)
      write(*, *) 'pv_pole_stokes      = ', to_str(pv_pole_stokes)
    if (pv_scheme == 'upwind') then
      write(*, *) 'upwind_order_pv     = ', to_str(upwind_order_pv)
      write(*, *) 'upwind_wgt_pv       = ', to_str(upwind_wgt_pv, 2)
    end if
      write(*, *) 'time_scheme         = ', trim(time_scheme)
      write(*, *) 'upwind_order        = ', to_str(upwind_order)
      write(*, *) 'use_topo_smooth     = ', to_str(use_topo_smooth)
    if (use_topo_smooth) then
      write(*, *) 'topo_max_slope      = ', topo_max_slope
      write(*, *) 'topo_smooth_cycles  = ', to_str(topo_smooth_cycles)
    end if
      write(*, *) 'use_div_damp        = ', to_str(use_div_damp)
    if (use_div_damp) then
      write(*, *) 'div_damp_cycles     = ', to_str(div_damp_cycles)
      write(*, *) 'div_damp_order      = ', to_str(div_damp_order)
      write(*, *) 'div_damp_coef2      = ', div_damp_coef2
      write(*, *) 'div_damp_coef4      = ', div_damp_coef4
      write(*, *) 'div_damp_top        = ', to_str(div_damp_top, 3)
      write(*, *) 'div_damp_pole_x     = ', to_str(div_damp_pole_x, 3)
      write(*, *) 'div_damp_pole_y     = ', to_str(div_damp_pole_y, 3)
      write(*, *) 'div_damp_lat0       = ', to_str(div_damp_lat0, 3)
    end if
    if (use_vor_damp) then
      write(*, *) 'vor_damp_cycles     = ', to_str(vor_damp_cycles)
      write(*, *) 'vor_damp_order      = ', to_str(vor_damp_order)
      write(*, *) 'vor_damp_coef2      = ', vor_damp_coef2
      write(*, *) 'vor_damp_k0         = ', to_str(vor_damp_k0)
      write(*, *) 'vor_damp_top        = ', to_str(vor_damp_top, 3)
      write(*, *) 'vor_damp_pole_x     = ', to_str(vor_damp_pole_x, 3)
      write(*, *) 'vor_damp_pole_y     = ', to_str(vor_damp_pole_y, 3)
      write(*, *) 'vor_damp_lat0       = ', to_str(vor_damp_lat0, 3)
    end if
    if (nonhydrostatic) then
      write(*, *) 'implicit_w_wgt      = ', to_str(implicit_w_wgt, 3)
      write(*, *) 'rayleigh_damp_w_coef= ', to_str(rayleigh_damp_w_coef, 2)
      write(*, *) 'rayleigh_damp_top   = ', to_str(rayleigh_damp_top   , 2)
    end if
      write(*, *) 'use_smag_damp       = ', to_str(use_smag_damp)
    if (use_smag_damp) then
      write(*, *) 'smag_damp_cycles    = ', to_str(smag_damp_cycles)
      write(*, *) 'smag_damp_coef      = ', smag_damp_coef
    end if
      write(*, *) '========================================================='

  end subroutine print_namelist

end module namelist_mod
