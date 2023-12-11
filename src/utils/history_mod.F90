! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module history_mod

  use mpi
  use container
  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod, dt => dt_dyn
  use time_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use allocator_mod
  use block_mod
  use tracer_mod
  use operators_mod, only: calc_div

  implicit none

  private

  public history_init_stage1
  public history_init_stage2
  public history_final
  public history_setup_h0_adv
  public history_write_h0
  public history_write_h1

  character(8), parameter ::    cell_dims_2d(3) = ['lon ', 'lat ',         'time']
  character(8), parameter ::    cell_dims_3d(4) = ['lon ', 'lat ', 'lev ', 'time']
  character(8), parameter ::     vtx_dims_2d(3) = ['ilon', 'ilat',         'time']
  character(8), parameter ::     vtx_dims_3d(4) = ['ilon', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lon_dims_2d(3) = ['ilon', 'lat ',         'time']
  character(8), parameter ::     lon_dims_3d(4) = ['ilon', 'lat ', 'lev ', 'time']
  character(8), parameter ::     lat_dims_2d(3) = ['lon ', 'ilat',         'time']
  character(8), parameter ::     lat_dims_3d(4) = ['lon ', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lev_dims_3d(4) = ['lon ', 'lat ', 'ilev', 'time']
  character(8), parameter :: lon_lev_dims_3d(4) = ['ilon', 'lat ', 'ilev', 'time']
  character(8), parameter :: lat_lev_dims_3d(4) = ['lon ', 'ilat', 'ilev', 'time']

contains

  subroutine history_init_stage1()

    character(10) time_value, time_units
    real(r8) seconds, months

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days', 'sol')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      if (proc%is_root()) call log_error('Invalid history interval ' // trim(history_interval(1)) // '!', __FILE__, __LINE__)
    end select

    call time_add_alert('history_write', seconds=seconds)

    if (trim(output_h0_new_file) == '') then
      call time_add_alert('h0_new_file', seconds=seconds)
      if (proc%is_root()) call log_notice('Output data every ' // trim(history_interval(1)) // '.')
    else if (output_h0_new_file == 'one_file') then
      if (proc%is_root()) call log_notice('Output data in one file.')
    else
      time_value = split_string(output_h0_new_file, ' ', 1)
      time_units = split_string(output_h0_new_file, ' ', 2)
      if (time_units == 'months') then
        read(time_value, *) months
        if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' months.')
        call time_add_alert('h0_new_file', months=months)
      else
        read(time_value, *) seconds
        select case (time_units)
        case ('days')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' days.')
          seconds = seconds * 86400
        case ('sol')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' sol.')
          seconds = seconds * 86400
        case ('hours')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' hours.')
          seconds = seconds * 3600
        case ('minutes')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' minutes.')
          seconds = seconds * 60
        case ('seconds')
          if (proc%is_root()) call log_notice('Output data every ' // trim(time_value) // ' seconds.')
          seconds = seconds
        case default
          if (proc%is_root()) call log_error('Invalid output_h0_new_file ' // trim(output_h0_new_file) // '!', __FILE__, __LINE__)
        end select
        call time_add_alert('h0_new_file', seconds=seconds)
      end if
    end if

    call fiona_init(time_units, start_time_str)

  end subroutine history_init_stage1

  subroutine history_init_stage2()

    if (hydrostatic) then
      call history_setup_h0_hydrostatic()
      call history_setup_h1_hydrostatic()
    else if (nonhydrostatic) then
      call history_setup_h0_nonhydrostatic()
      call history_setup_h1_nonhydrostatic()
    else if (.not. advection) then
      call history_setup_h0_swm()
      call history_setup_h1_swm()
    end if

  end subroutine history_init_stage2

  subroutine history_final()

    call fiona_final()

  end subroutine history_final

  subroutine history_setup_h0_swm()

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Global attributes
    call fiona_add_att('h0', 'planet', planet)
    call fiona_add_att('h0', 'time_step_size', dt)
    ! Dimensions
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    ! Variables
    select case (planet)
    case ('mars')
      call fiona_add_var('h0', 'Ls', long_name='Solar longitude', units='deg', dim_names=['time'], dtype=output_h0_dtype)
    end select
    call fiona_add_var('h0', 'tm'   , long_name='Total mass'               , units='m'      , dim_names=['time'])
    call fiona_add_var('h0', 'te'   , long_name='Total energy'             , units='m4 s-4' , dim_names=['time']    , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'tpe'  , long_name='Total potential enstrophy', units='m2 s-5' , dim_names=['time']    , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'    , long_name='U wind component'         , units='m s-1'  , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'    , long_name='V wind component'         , units='m s-1'  , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'z'    , long_name='Height'                   , units='m'      , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pv'   , long_name='Potential vorticity'      , units='s-1 m-1', dim_names= vtx_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'div'  , long_name='Divergence'               , units='s-1'    , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    if (test_case == 'stratospheric_vortex_erosion') then
      call fiona_add_var('h0', 'zs', long_name='Surface height', units='m', dim_names=cell_dims_2d)
    else
      call fiona_add_var('h0', 'zs', long_name='Surface height', units='m', dim_names=['lon', 'lat'])
    end if

  end subroutine history_setup_h0_swm

  subroutine history_setup_h0_adv(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer i, j

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h0', 'z', long_name='Geopotential height', units='m', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    do i = 1, ntracers
      call fiona_add_var('h0', tracer_names(i), long_name=tracer_long_names(i), units=tracer_units(i), dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end do

  end subroutine history_setup_h0_adv

  subroutine history_setup_h0_hydrostatic()

    integer k

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), &
      mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Global attributes
    call fiona_add_att('h0', 'planet', planet)
    call fiona_add_att('h0', 'time_step_size', dt)
    ! Dimensions
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
    ! Variables
    select case (planet)
    case ('mars')
      call fiona_add_var('h0', 'Ls', long_name='solar longitude', units='deg', dim_names=['time'], dtype=output_h0_dtype)
    end select
    call fiona_add_var('h0', 'tm'     , long_name='Total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='Total energy'                , units='m4 s-4', dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'tpe'    , long_name='Total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_ke'  , long_name='Total kinetic energy'        , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_ie'  , long_name='Total internal energy'       , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_pe'  , long_name='Total potential energy'      , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'zs'     , long_name='Surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    call fiona_add_var('h0', 'phs'    , long_name='Surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'      , long_name='U wind component'            , units='m s-1' , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'      , long_name='V wind component'            , units='m s-1' , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pt'     , long_name='Potential temperature'       , units='K'     , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 't'      , long_name='Temperature'                 , units='K'     , dim_names=cell_dims_3d  )
    call fiona_add_var('h0', 'z'      , long_name='Height'                      , units='m'     , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ph'     , long_name='Hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'vor'    , long_name='Relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'div'    , long_name='Divergence'                  , units='s-1'   , dim_names=cell_dims_3d  )
    if (use_div_damp .and. div_damp_order == 4) then
      call fiona_add_var('h0', 'div2' , long_name='Laplacian of divergence', units='s-1 m-2', dim_names=cell_dims_3d)
    end if
    call fiona_add_var('h0', 'landmask', long_name='Land mask'                  , units=''      , dim_names=['lon', 'lat'])

    if (vert_coord_scheme == 'smooth') then
      call fiona_add_var('h0', 'ref_ps', long_name='Reference surface pressure', units='Pa', dim_names=['lon', 'lat'], dtype=output_h0_dtype)
      call fiona_add_var('h0', 'ref_ps_smth', long_name='Smoothed reference surface pressure', units='Pa', dim_names=['lon', 'lat'], dtype=output_h0_dtype)
      call fiona_add_var('h0', 'ref_ps_perb', long_name='Reference surface pressure perturbation', units='Pa', dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    end if

    if (idx_qv > 0) then
      call fiona_add_var('h0', 'qv', long_name='Water vapor mixing ratio', units='kg kg-1', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end if
    if (idx_qc > 0) then
      call fiona_add_var('h0', 'qc', long_name='Cloud water mixing ratio', units='kg kg-1', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end if
    if (idx_nc > 0) then
      call fiona_add_var('h0', 'nc', long_name='Cloud water number concentration', units='m-3', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end if
    if (idx_qi > 0) then
      call fiona_add_var('h0', 'qi', long_name='Cloud ice mixing ratio', units='kg kg-1', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end if
    if (idx_ni > 0) then
      call fiona_add_var('h0', 'ni', long_name='Cloud ice number concentration', units='m-3', dim_names=cell_dims_3d, dtype=output_h0_dtype)
    end if

    do k = 1, blocks(1)%accum_list%size
      select type (accum => blocks(1)%accum_list%value_at(k))
      type is (accum_type)
        if (.not. accum%active) cycle
        call fiona_add_var('h0', accum%name, accum%units, accum%long_name, dim_names=cell_dims_3d, dtype=output_h0_dtype)
      end select
    end do

    if (test_case == 'tropical_cyclone') then
      call fiona_add_var('h0', 'precl', '', 'Large-scale precipitation', dim_names=cell_dims_2d, dtype=output_h0_dtype)
    end if

  end subroutine history_setup_h0_hydrostatic

  subroutine history_setup_h0_nonhydrostatic()

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h0', 'tm'     , long_name='Total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='Total energy'                , units='m4 s-4', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'tpe'    , long_name='Total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'zs'     , long_name='Surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype='r8')
    call fiona_add_var('h0', 'phs'    , long_name='Surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'      , long_name='U wind component'            , units='m s-1' , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'      , long_name='V wind component'            , units='m s-1' , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pt'     , long_name='Potential temperature'       , units='K'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 't'      , long_name='Temperature'                 , units='K'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'z'      , long_name='Height'                      , units='m'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'z_lev'  , long_name='Height'                      , units='m'     , dim_names= lev_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ph'     , long_name='Hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'vor'    , long_name='Relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d)
    call fiona_add_var('h0', 'div'    , long_name='Divergence'                  , units='s-1'   , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'w_lev'  , long_name='Vertical speed'              , units='m s-1' , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'p'      , long_name='Pressure'                    , units='Pa'    , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'rhod'   , long_name='Dry-air density'             , units='kg m-3', dim_names=cell_dims_3d)

  end subroutine history_setup_h0_nonhydrostatic

  subroutine history_setup_h1_swm()

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    ! Variables
    call fiona_add_var('h1', 'dudt'    , long_name='U wind component tendency'                     , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'dvdt'    , long_name='V wind component tendency'                     , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'dgzdt'   , long_name='Geopotential tendency'                         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dmf'     , long_name='Mass flux divergence'                          , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'mfx_lon' , long_name='Normal mass flux on U grid'                    , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'mfy_lat' , long_name='Normal mass flux on V grid'                    , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'ke'      , long_name='Kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d)

  end subroutine history_setup_h1_swm

  subroutine history_setup_h1_hydrostatic()

    integer m
    character(30) tag

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lev'  , size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilev' , size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h1', 'u_lon'        , long_name='U wind component'                              , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'v_lat'        , long_name='V wind component'                              , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dudt'         , long_name='U wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='V wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmgsdt'       , long_name='Surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='Potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev'       , long_name='Vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'gz_lev'       , long_name='Geopotential on half levels'                   , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'ph_lev'       , long_name='Hydrostatic pressure on half levels'           , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'dmf'          , long_name='Mass flux divergence'                          , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'omg'          , long_name='Vertical pressure velocity'                    , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='Normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='Normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmg'          , long_name='Dry-air weight on full levels'                 , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dmg_lon'      , long_name='Dry-air weight on U grid'                      , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dmg_lat'      , long_name='Dry-air weight on V grid'                      , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmg_vtx'      , long_name='Dry-air weight on PV grid'                     , units='', dim_names= vtx_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='Kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'tv'           , long_name='Virtual temperature'                           , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'n2_lev'       , long_name='Square of buoyancy frequency'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'ri_lev'       , long_name='Local Richardson number'                       , units='', dim_names= lev_dims_3d)

#ifdef OUTPUT_H1_DTEND
    call fiona_add_var('h1', 'dudt_coriolis', long_name='Nonlinear Coriolis tendency'                   , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt_coriolis', long_name='Nonlinear Coriolis tendency'                   , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dudt_wedudeta', long_name='Vertical advection tendency of U'              , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt_wedvdeta', long_name='Vertical advection tendency of V'              , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dudt_dkedx'   , long_name='Kinetic energy gradient tendency of U'         , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt_dkedy'   , long_name='Kinetic energy gradient tendency of V'         , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dudt_pgf'     , long_name='Horizontal PGF tendency of U'                  , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt_pgf'     , long_name='Horizontal PGF tendency of V'                  , units='', dim_names= lat_dims_3d)
#endif

    if (physics_suite /= 'N/A') then
      call fiona_add_var('h1', 'dudt_phys'  , long_name='Physics tendency for u'                        , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dvdt_phys'  , long_name='Physics tendency for v'                        , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dptdt_phys' , long_name='Physics tendency for pt'                       , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dqdt_phys'  , long_name='Physics tendency for q'                        , units='', dim_names=cell_dims_3d)
    end if

    do m = 1, size(blocks(1)%adv_batches)
      tag = blocks(1)%adv_batches(m)%name
      call fiona_add_var('h1', trim(tag)//'_mfx' , long_name='', units='', dim_names=lon_dims_3d)
      call fiona_add_var('h1', trim(tag)//'_mfy' , long_name='', units='', dim_names=lat_dims_3d)
      call fiona_add_var('h1', trim(tag)//'_cflx', long_name='', units='', dim_names=lon_dims_3d)
      call fiona_add_var('h1', trim(tag)//'_cfly', long_name='', units='', dim_names=lat_dims_3d)
      call fiona_add_var('h1', trim(tag)//'_cflz', long_name='', units='', dim_names=lev_dims_3d)
    end do

  end subroutine history_setup_h1_hydrostatic

  subroutine history_setup_h1_nonhydrostatic

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, ngroup=output_ngroup)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lev'  , size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilev' , size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h1', 'dudt'         , long_name='U wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='V wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmgsdt'       , long_name='Surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='Potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev'       , long_name='Vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'dmf'          , long_name='Divergence flux'                               , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='Normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='Normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmg'          , long_name='Dry-air weight on full levels'                 , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='Kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'adv_w_lev'    , long_name='Advection of w'                                , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_gz_lev'   , long_name='Advection of gz'                               , units='', dim_names= lev_dims_3d)

  end subroutine history_setup_h1_nonhydrostatic

  subroutine history_write_h0_swm(itime)

    integer, intent(in) :: itime

    integer iblk, is, ie, js, je
    integer start(2), count(2)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 aux    => blocks(iblk)%aux          , &
                 static => blocks(iblk)%static)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      start = [is,js]
      count = [mesh%full_nlon,mesh%full_nlat]
      call fiona_output('h0', 'u'  , dstate%u  %d(is:ie,js:je,1)    , start=start, count=count)
      call fiona_output('h0', 'v'  , dstate%v  %d(is:ie,js:je,1)    , start=start, count=count)
      call fiona_output('h0', 'zs' , static%gzs%d(is:ie,js:je  ) / g, start=start, count=count)
      call fiona_output('h0', 'z'  , dstate%gz %d(is:ie,js:je,1) / g, start=start, count=count)
      call fiona_output('h0', 'div', aux%div   %d(is:ie,js:je,1)    , start=start, count=count)
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      start = [is,js]
      count = [mesh%half_nlon,mesh%half_nlat]
      call fiona_output('h0', 'pv' , aux%pv    %d(is:ie,js:je,1)    , start=start, count=count)

      call fiona_output('h0', 'tm' , dstate %tm)
      call fiona_output('h0', 'te' , dstate %te)
      call fiona_output('h0', 'tpe', dstate %tpe)
      end associate
    end do

  end subroutine history_write_h0_swm

  subroutine history_write_h0_adv(itime)

    integer, intent(in) :: itime

    integer iblk, is, ie, js, je, ks, ke, i
    integer start(3), count(3)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h0', 'z', blocks(1)%dstate(itime)%gz%d(is:ie,js:je,ks:ke) / g, start=start, count=count)
      do i = 1, ntracers
        call fiona_output('h0', tracer_names(i), tracers(iblk)%q%d(is:ie,js:je,ks:ke,i), start=start, count=count)
      end do
      end associate
    end do

  end subroutine history_write_h0_adv

  subroutine history_write_h0_hydrostatic(itime)

    integer, intent(in) :: itime

    integer iblk, is, ie, js, je, ks, ke, k
    integer start(3), count(3)
    real(r8), pointer :: q(:,:,:)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))

    do iblk = 1, size(blocks)
      associate (mesh        => blocks(iblk)%mesh             , &
                 gzs         => blocks(iblk)%static%gzs       , &
                 landmask    => blocks(iblk)%static%landmask  , &
                 tm          => blocks(iblk)%dstate(itime)%tm , &
                 te          => blocks(iblk)%dstate(itime)%te , &
                 u           => blocks(iblk)%dstate(itime)%u  , &
                 v           => blocks(iblk)%dstate(itime)%v  , &
                 gz          => blocks(iblk)%dstate(itime)%gz , &
                 phs         => blocks(iblk)%dstate(itime)%phs, &
                 ph          => blocks(iblk)%dstate(itime)%ph , &
                 pt          => blocks(iblk)%dstate(itime)%pt , &
                 t           => blocks(iblk)%dstate(itime)%t  , &
                 div         => blocks(iblk)%aux%div          , &
                 div2        => blocks(iblk)%aux%div2         , &
                 vor         => blocks(iblk)%aux%vor          , &
                 q           => tracers(iblk)%q               , &
                 precl       => blocks(iblk)%pstate%precl     , &
                 accum_list  => blocks(iblk)%accum_list       )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h0', 'zs'      , gzs     %d(is:ie,js:je) / g  , start=start, count=count)
      call fiona_output('h0', 'landmask', landmask%d(is:ie,js:je)      , start=start, count=count)
      call fiona_output('h0', 'u'       , u       %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h0', 'v'       , v       %d(is:ie,js:je,ks:ke), start=start, count=count)
      gz%d = gz%d / g
      call fiona_output('h0', 'z'       , gz      %d(is:ie,js:je,ks:ke), start=start, count=count)
      gz%d = gz%d * g
      call fiona_output('h0', 'phs'     , phs     %d(is:ie,js:je)      , start=start, count=count)
      call fiona_output('h0', 'ph'      , ph      %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h0', 'pt'      , pt      %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h0', 't'       , t       %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h0', 'div'     , div     %d(is:ie,js:je,ks:ke), start=start, count=count)
      if (use_div_damp .and. div_damp_order == 4) then
        call fiona_output('h0', 'div2'  , div2    %d(is:ie,js:je,ks:ke) , start=start, count=count)
      end if
      if (idx_qv > 0) then
        call fiona_output('h0', 'qv'    , q%d(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
      end if
      if (idx_qc > 0) then
        call fiona_output('h0', 'qc'    , q%d(is:ie,js:je,ks:ke,idx_qc), start=start, count=count)
      end if
      if (idx_nc > 0) then
        call fiona_output('h0', 'nc'    , q%d(is:ie,js:je,ks:ke,idx_nc), start=start, count=count)
      end if
      if (idx_qi > 0) then
        call fiona_output('h0', 'qi'    , q%d(is:ie,js:je,ks:ke,idx_qi), start=start, count=count)
      end if
      if (idx_ni > 0) then
        call fiona_output('h0', 'ni'    , q%d(is:ie,js:je,ks:ke,idx_ni), start=start, count=count)
      end if
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_output('h0', 'vor', vor%d(is:ie,js:je,ks:ke), start=start, count=count)

      call fiona_output('h0', 'tm', tm)
      call fiona_output('h0', 'te', te)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      do k = 1, accum_list%size
        select type (accum => accum_list%value_at(k))
        type is (accum_type)
          if (.not. accum%active) cycle
          call fiona_output('h0', accum%name, accum%array(:,:,:,1), start=start, count=count)
        end select
      end do

      if (test_case == 'tropical_cyclone') then
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        call fiona_output('h0', 'precl', reshape(precl, count(1:2)), start=start(1:2), count=count(1:2))
      end if
      end associate
    end do

  end subroutine history_write_h0_hydrostatic

  subroutine history_write_h0_nonhydrostatic(itime)

    integer, intent(in) :: itime

    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 aux    => blocks(iblk)%aux          , &
                 static => blocks(iblk)%static)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h0', 'zs'      , static%gzs    %d(is:ie,js:je) / g       , start=start, count=count)
      call fiona_output('h0', 'u'       , dstate%u      %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'v'       , dstate%v      %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'z'       , dstate%gz     %d(is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'phs'     , dstate%phs    %d(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'ph'      , dstate%ph     %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'pt'      , dstate%pt     %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 't'       , dstate%t      %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'div'     , aux%div       %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'rhod'    , dstate%rhod   %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_output('h0', 'vor'     , aux%vor       %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
      call fiona_output('h0', 'w_lev'   , dstate%w_lev   %d(is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'z_lev'   , dstate%gz_lev  %d(is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'p'       , dstate%p_lev   %d(is:ie,js:je,ks:ke)     , start=start, count=count)

      call fiona_output('h0', 'tm' , dstate %tm)
      call fiona_output('h0', 'te' , dstate %te)
      call fiona_output('h0', 'tpe', dstate %tpe)
      end associate
    end do

  end subroutine history_write_h0_nonhydrostatic

  subroutine history_write_h0(itime)

    integer, intent(in) :: itime

    logical, save :: first_call = .true.
    real(8) time1, time2

    if (proc%is_root()) then
      call log_notice('Write history.')
      time1 = MPI_WTIME()
    end if

    if (.not. use_div_damp) then
      call calc_div(blocks(1), blocks(1)%dstate(itime))
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=first_call)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=first_call)
    end if

    first_call = .false.

    select case (planet)
    case ('mars')
      call fiona_output('h0', 'Ls', curr_time%solar_longitude() * deg)
    end select

    if (advection) then
      call history_write_h0_adv(itime)
    else if (hydrostatic) then
      call history_write_h0_hydrostatic(itime)
    else if (nonhydrostatic) then
      call history_write_h0_nonhydrostatic(itime)
    else
      call history_write_h0_swm(itime)
    end if

    call fiona_end_output('h0', keep_dataset=.true.)

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write history cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_h0

  subroutine history_write_h1_swm(itime)

    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke
    integer start(2), count(2)

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%half_nlat))

    associate (mesh   => blocks(1)%mesh         , &
               dstate => blocks(1)%dstate(itime), &
               dtend  => blocks(1)%dtend (itime), &
               aux    => blocks(1)%aux          , &
               static => blocks(1)%static)
    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    start = [is,js]
    count = [mesh%full_nlon,mesh%full_nlat]
    call fiona_output('h1', 'dmf'     ,    aux%dmf     %d(is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke      %d(is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dgzdt'   ,  dtend%dgz     %d(is:ie,js:je,1), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    start = [is,js]
    count = [mesh%half_nlon,mesh%full_nlat]
    call fiona_output('h1', 'dudt   ' ,  dtend%du      %d(is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' ,    aux%mfx_lon %d(is:ie,js:je,1), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    start = [is,js]
    count = [mesh%full_nlon,mesh%half_nlat]
    call fiona_output('h1', 'dvdt'    ,  dtend%dv      %d(is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' ,    aux%mfy_lat %d(is:ie,js:je,1), start=start, count=count)
    end associate

    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_swm

  subroutine history_write_h1_hydrostatic(itime)

    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke, m
    integer start(3), count(3)
    character(30) tag

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h1', 'lev'   , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h1', 'ilev'  , global_mesh%half_lev(1:global_mesh%half_nlev))

    associate (mesh   => blocks(1)%mesh         , &
               dstate => blocks(1)%dstate(itime), &
               pstate => blocks(1)%pstate       , &
               dtend  => blocks(1)%dtend (itime), &
               aux    => blocks(1)%aux          , &
               static => blocks(1)%static)
    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dmf'     ,    aux%dmf    %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'omg'     ,    aux%omg    %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke     %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'tv'      , dstate%tv     %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg'     , dstate%dmg    %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmgsdt'  ,  dtend%dmgs   %d(is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  dtend%dpt    %d(is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dmg_vtx' ,    aux%dmg_vtx%d(is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'u_lon'   , dstate%u_lon  %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt   ' ,  dtend%du     %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' ,    aux%mfx_lon%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg_lon' ,    aux%dmg_lon%d(is:ie,js:je,ks:ke), start=start, count=count)
#ifdef OUTPUT_H1_DTEND
    call fiona_output('h1', 'dudt_coriolis', dtend%dudt_coriolis%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt_wedudeta', dtend%dudt_wedudeta%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt_dkedx'   , dtend%dudt_dkedx   %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt_pgf'     , dtend%dudt_pgf     %d(is:ie,js:je,ks:ke), start=start, count=count)
#endif
    do m = 1, size(blocks(1)%adv_batches)
      tag = blocks(1)%adv_batches(m)%name
      call fiona_output('h1', trim(tag)//'_mfx' , blocks(1)%adv_batches(m)%mfx %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', trim(tag)//'_cflx', blocks(1)%adv_batches(m)%cflx%d(is:ie,js:je,ks:ke), start=start, count=count)
    end do

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
    call fiona_output('h1', 'v_lat'   , dstate%v_lat    %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt'    ,  dtend%dv       %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' ,    aux%mfy_lat  %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg_lat' ,    aux%dmg_lat  %d(is:ie,js:je,ks:ke), start=start, count=count)
#ifdef OUTPUT_H1_DTEND
    call fiona_output('h1', 'dvdt_coriolis', dtend%dvdt_coriolis%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt_wedvdeta', dtend%dvdt_wedvdeta%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt_dkedy'   , dtend%dvdt_dkedy   %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt_pgf'     , dtend%dvdt_pgf     %d(is:ie,js:je,ks:ke), start=start, count=count)
#endif
    do m = 1, size(blocks(1)%adv_batches)
      tag = blocks(1)%adv_batches(m)%name
      call fiona_output('h1', trim(tag)//'_mfy' , blocks(1)%adv_batches(m)%mfy %d(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', trim(tag)//'_cfly', blocks(1)%adv_batches(m)%cfly%d(is:ie,js:je,ks:ke), start=start, count=count)
    end do

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%half_kds; ke = mesh%half_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
    call fiona_output('h1', 'we_lev', dstate%we_lev%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'gz_lev', dstate%gz_lev%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ph_lev', dstate%ph_lev%d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'n2_lev', reshape(pstate%n2_lev, count), start=start, count=count)
    call fiona_output('h1', 'ri_lev', reshape(pstate%ri_lev, count), start=start, count=count)
    do m = 1, size(blocks(1)%adv_batches)
      tag = blocks(1)%adv_batches(m)%name
      call fiona_output('h1', trim(tag)//'_cflz', blocks(1)%adv_batches(m)%cflz%d(is:ie,js:je,ks:ke), start=start, count=count)
    end do

    if (physics_suite /= 'N/A') then
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h1',  'dudt_phys', aux% dudt_phys%d(is:ie,js:je,ks:ke       ), start=start, count=count)
      call fiona_output('h1',  'dvdt_phys', aux% dvdt_phys%d(is:ie,js:je,ks:ke       ), start=start, count=count)
      call fiona_output('h1', 'dptdt_phys', aux%dptdt_phys%d(is:ie,js:je,ks:ke       ), start=start, count=count)
      call fiona_output('h1',  'dqdt_phys', aux% dqdt_phys%d(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
    end if
    end associate

    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_hydrostatic

  subroutine history_write_h1_nonhydrostatic(itime)

    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', dble(elapsed_seconds), new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h1', 'lev'   , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h1', 'ilev'  , global_mesh%half_lev(1:global_mesh%half_nlev))

    associate (mesh   => blocks(1)%mesh         , &
               dstate => blocks(1)%dstate(itime), &
               dtend  => blocks(1)%dtend (itime), &
               aux    => blocks(1)%aux          , &
               static => blocks(1)%static)
    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dmf'     ,    aux%dmf      %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke       %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg'     , dstate%dmg      %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmgsdt'  ,  dtend%dmgs     %d(is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  dtend%dpt      %d(is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dudt   ' ,  dtend%du       %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' ,    aux%mfx_lon  %d(is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dvdt'    ,  dtend%dv       %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' ,    aux%mfy_lat  %d(is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%half_kds; ke = mesh%half_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
    call fiona_output('h1', 'we_lev'    , dstate%we_lev    %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_w_lev' ,    aux%adv_w_lev %d(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_gz_lev',    aux%adv_gz_lev%d(is:ie,js:je,ks:ke), start=start, count=count)
    end associate

    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_nonhydrostatic

  subroutine history_write_h1(itime)

    integer, intent(in) :: itime

    if (hydrostatic) then
      call history_write_h1_hydrostatic(itime)
    else if (nonhydrostatic) then
      call history_write_h1_nonhydrostatic(itime)
    else
      call history_write_h1_swm(itime)
    end if

  end subroutine history_write_h1

end module history_mod
