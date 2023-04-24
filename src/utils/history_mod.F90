module history_mod

  use mpi
  use container
  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod, dt => dt_dyn
  use time_mod
  use parallel_mod
  use allocator_mod
  use block_mod
  use diag_state_mod
  use tracer_mod

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
      call fiona_add_var('h0', 'Ls', long_name='solar longitude', units='deg', dim_names=['time'], dtype=output_h0_dtype)
    end select
    call fiona_add_var('h0', 'tm'   , long_name='total mass'               , units='m'      , dim_names=['time'])
    call fiona_add_var('h0', 'te'   , long_name='total energy'             , units='m4 s-4' , dim_names=['time']    , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy', units='m2 s-5' , dim_names=['time']    , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'    , long_name='u wind component'         , units='m s-1'  , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'    , long_name='v wind component'         , units='m s-1'  , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'z'    , long_name='height'                   , units='m'      , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'      , units='s-1 m-1', dim_names= vtx_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'div'  , long_name='divergence'               , units='s-1'    , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    if (test_case == 'stratospheric_vortex_erosion') then
      call fiona_add_var('h0', 'zs', long_name='surface height', units='m', dim_names=cell_dims_2d)
    else
      call fiona_add_var('h0', 'zs', long_name='surface height', units='m', dim_names=['lon', 'lat'])
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
    call fiona_add_var('h0', 'tm'     , long_name='total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='total energy'                , units='m4 s-4', dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'tpe'    , long_name='total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_ke'  , long_name='total kinetic energy'        , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_ie'  , long_name='total internal energy'       , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'te_pe'  , long_name='total potential energy'      , units=''      , dim_names=['time']      , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ref_ps' , long_name='reference surface pressure'  , units='Pa'    , dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ref_ps_smth' , long_name='smoothed reference surface pressure'  , units='Pa'    , dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ref_ps_perb' , long_name='reference surface pressure perturbation'  , units='Pa'    , dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    call fiona_add_var('h0', 'zs'     , long_name='surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype=output_h0_dtype)
    call fiona_add_var('h0', 'dzsdlon', long_name='zonal zs gradient'           , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'dzsdlat', long_name='meridional zs gradient'      , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'phs'    , long_name='surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'      , long_name='u wind component'            , units='m s-1' , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'      , long_name='v wind component'            , units='m s-1' , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pt'     , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 't'      , long_name='temperature'                 , units='K'     , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'z'      , long_name='height'                      , units='m'     , dim_names=cell_dims_3d  , dtype=output_h0_dtype)
    ! call fiona_add_var('h0', 'ph'     , long_name='hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'vor'    , long_name='relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d)
    ! call fiona_add_var('h0', 'div'    , long_name='divergence'                  , units='s-1'   , dim_names=cell_dims_3d)
    ! call fiona_add_var('h0', 'landmask', long_name='land mask'                  , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'qv'     , long_name='water vapor mixing ratio'    , units='g g-1' , dim_names=cell_dims_3d  , dtype=output_h0_dtype)

    select case (diag_state(1)%level_type)
    case (height_levels)
      do k = 1, size(diag_state(1)%levels)
        call fiona_add_var('h0', 'u' // to_str(int(diag_state(1)%levels(k))) // 'm', &
          long_name='zonal wind speed on ' // to_str(int(diag_state(1)%levels(k))) // 'm', &
          units='m s-1', dim_names=cell_dims_2d)
        call fiona_add_var('h0', 'v' // to_str(int(diag_state(1)%levels(k))) // 'm', &
          long_name='meridional wind speed on ' // to_str(int(diag_state(1)%levels(k))) // 'm', &
          units='m s-1', dim_names=cell_dims_2d)
      end do
    case (pressure_levels)
      do k = 1, size(diag_state(1)%levels)
        call fiona_add_var('h0', 'u' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          long_name='zonal wind speed on ' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          units='m s-1', dim_names=cell_dims_2d)
        call fiona_add_var('h0', 'v' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          long_name='meridional wind speed on ' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          units='m s-1', dim_names=cell_dims_2d)
        call fiona_add_var('h0', 't' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          long_name='temperature on ' // to_str(int(diag_state(1)%levels(k))/100) // 'hPa', &
          units='K', dim_names=cell_dims_2d)
      end do
    end select

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
    call fiona_add_var('h0', 'tm'     , long_name='total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='total energy'                , units='m4 s-4', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'tpe'    , long_name='total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'zs'     , long_name='surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype='r8')
    call fiona_add_var('h0', 'dzsdlon', long_name='zonal zs gradient'           , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'dzsdlat', long_name='meridional zs gradient'      , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'phs'    , long_name='surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'u'      , long_name='u wind component'            , units='m s-1' , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'v'      , long_name='v wind component'            , units='m s-1' , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'pt'     , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 't'      , long_name='temperature'                 , units='K'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'z'      , long_name='height'                      , units='m'     , dim_names=cell_dims_3d, dtype=output_h0_dtype)
    call fiona_add_var('h0', 'ph'     , long_name='hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'vor'    , long_name='relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d)
    call fiona_add_var('h0', 'div'    , long_name='divergence'                  , units='s-1'   , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'w'      , long_name='vertical speed'              , units='m s-1' , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'p'      , long_name='full pressure'               , units='Pa'    , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'rhod'   , long_name='dry air density'             , units='kg m-3', dim_names=cell_dims_3d)

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
    call fiona_add_var('h1', 'dudt'    , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'dvdt'    , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'dgzdt'   , long_name='geopotential tendency'                         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dmf'     , long_name='mass flux divergence'                          , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'mfx_lon' , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'mfy_lat' , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'ke'      , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d)

  end subroutine history_setup_h1_swm

  subroutine history_setup_h1_hydrostatic()

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
    call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmgsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev'       , long_name='vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'dmf'          , long_name='mass flux divergence'                          , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmg'          , long_name='dry-air weight on full levels'                 , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'n2'           , long_name='square of buoyancy frequency'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'ri'           , long_name='local Richardson number'                       , units='', dim_names= lev_dims_3d)

    if (physics_suite /= 'none') then
      call fiona_add_var('h1', 'dudt_phys'  , long_name='physics tendency for u'                        , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dvdt_phys'  , long_name='physics tendency for v'                        , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dtdt_phys'  , long_name='physics tendency for t'                        , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'dshdt_phys' , long_name='physics tendency for sh'                       , units='', dim_names=cell_dims_3d)
    end if

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
    call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmgsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev'       , long_name='vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'dmf'          , long_name='flux divergence'                               , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmg'          , long_name='dry-air weight on full levels'                 , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'adv_gz'       , long_name='advection of gz'                               , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_w'        , long_name='advection of w'                                , units='', dim_names= lev_dims_3d)

  end subroutine history_setup_h1_nonhydrostatic

  subroutine history_write_h0_swm(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
      call fiona_output('h0', 'u'  , dstate%u  (is:ie,js:je,1)    , start=start, count=count)
      call fiona_output('h0', 'v'  , dstate%v  (is:ie,js:je,1)    , start=start, count=count)
      call fiona_output('h0', 'zs' , static%gzs(is:ie,js:je  ) / g, start=start, count=count)
      call fiona_output('h0', 'z'  , dstate%gz (is:ie,js:je,1) / g, start=start, count=count)
      call fiona_output('h0', 'div', aux%div   (is:ie,js:je,1)    , start=start, count=count)
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      start = [is,js]
      count = [mesh%half_nlon,mesh%half_nlat]
      call fiona_output('h0', 'pv' , aux%pv    (is:ie,js:je,1)    , start=start, count=count)
 
      call fiona_output('h0', 'tm' , dstate %tm)
      call fiona_output('h0', 'te' , dstate %te)
      call fiona_output('h0', 'tpe', dstate %tpe)
      end associate
    end do

  end subroutine history_write_h0_swm

  subroutine history_write_h0_adv(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
      do i = 1, ntracers
        call fiona_output('h0', tracer_names(i), tracers(iblk)%q(is:ie,js:je,ks:ke,i), start=start, count=count)
      end do
      end associate
    end do

  end subroutine history_write_h0_adv

  subroutine history_write_h0_hydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime 

    integer iblk, is, ie, js, je, ks, ke, k
    integer start(3), count(3)
    real(r8), pointer :: qv(:,:,:)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%full_nlev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))

    do iblk = 1, size(blocks)
      associate (mesh        => blocks(iblk)%mesh         , &
                 dstate      => blocks(iblk)%dstate(itime), &
                 static      => blocks(iblk)%static       , &
                 adv_batches => blocks(iblk)%adv_batches  , &
                 aux         => blocks(iblk)%aux          , &
                 pstate      => blocks(iblk)%pstate       )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h0', 'ref_ps'  , static%ref_ps (is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'ref_ps_smth', static%ref_ps_smth(is:ie,js:je)    , start=start, count=count)
      call fiona_output('h0', 'ref_ps_perb', static%ref_ps_perb(is:ie,js:je)    , start=start, count=count)
      call fiona_output('h0', 'zs'      , static%gzs    (is:ie,js:je) / g       , start=start, count=count)
      call fiona_output('h0', 'dzsdlon' , static%dzsdlon(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'dzsdlat' , static%dzsdlat(is:ie,js:je)           , start=start, count=count)
      ! call fiona_output('h0', 'landmask', static%landmask(is:ie,js:je)          , start=start, count=count)
      call fiona_output('h0', 'u'       , dstate%u      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'v'       , dstate%v      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'z'       , dstate%gz     (is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'phs'     , dstate%phs    (is:ie,js:je)           , start=start, count=count)
      ! call fiona_output('h0', 'ph'      , dstate%ph     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'pt'      , dstate%pt     (is:ie,js:je,ks:ke)     , start=start, count=count)
      ! call fiona_output('h0', 't'       , dstate%t      (is:ie,js:je,ks:ke)     , start=start, count=count)
      ! call fiona_output('h0', 'div'     , aux%div       (is:ie,js:je,ks:ke)     , start=start, count=count)
      call tracer_get_array(iblk, idx_qv, qv)
      call fiona_output('h0', 'qv'      , qv            (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_output('h0', 'vor'     , aux%vor       (is:ie,js:je,ks:ke)     , start=start, count=count)
      
      call fiona_output('h0', 'tm'   , dstate%tm)
      call fiona_output('h0', 'te'   , dstate%te)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      select case (diag_state(1)%level_type)
      case (height_levels)
        do k = 1, size(diag_state(iblk)%levels)
          call fiona_output('h0', 'u' // to_str(int(diag_state(iblk)%levels(k))) // 'm', &
            diag_state(iblk)%u(is:ie,js:je,k), start=start, count=count)
          call fiona_output('h0', 'v' // to_str(int(diag_state(iblk)%levels(k))) // 'm', &
            diag_state(iblk)%v(is:ie,js:je,k), start=start, count=count)
        end do
      case (pressure_levels)
        do k = 1, size(diag_state(iblk)%levels)
          call fiona_output('h0', 'u' // to_str(int(diag_state(iblk)%levels(k))/100) // 'hPa', &
            diag_state(iblk)%u(is:ie,js:je,k), start=start, count=count)
          call fiona_output('h0', 'v' // to_str(int(diag_state(iblk)%levels(k))/100) // 'hPa', &
            diag_state(iblk)%v(is:ie,js:je,k), start=start, count=count)
          call fiona_output('h0', 't' // to_str(int(diag_state(iblk)%levels(k))/100) // 'hPa', &
            diag_state(iblk)%t(is:ie,js:je,k), start=start, count=count)
        end do
      end select
      end associate
    end do

  end subroutine history_write_h0_hydrostatic

  subroutine history_write_h0_nonhydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
      call fiona_output('h0', 'zs'      , static%gzs    (is:ie,js:je) / g       , start=start, count=count)
      call fiona_output('h0', 'dzsdlon' , static%dzsdlon(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'dzsdlat' , static%dzsdlat(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'u'       , dstate%u      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'v'       , dstate%v      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'z'       , dstate%gz     (is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'phs'     , dstate%phs    (is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'ph'      , dstate%ph     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'pt'      , dstate%pt     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 't'       , dstate%t      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'div'     , aux%div       (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'rhod'    , dstate%rhod   (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_output('h0', 'vor'     , aux%vor       (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
      call fiona_output('h0', 'w'       , dstate%w_lev   (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'p'       , dstate%p_lev   (is:ie,js:je,ks:ke)     , start=start, count=count)
      
      call fiona_output('h0', 'tm' , dstate %tm)
      call fiona_output('h0', 'te' , dstate %te)
      call fiona_output('h0', 'tpe', dstate %tpe)
      end associate
    end do

  end subroutine history_write_h0_nonhydrostatic

  subroutine history_write_h0(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    real(8) time1, time2

    if (proc%is_root()) then
      call log_notice('Write state.')
      time1 = MPI_WTIME()
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h0', dble(elapsed_seconds), new_file=time_step==0)
    end if

    select case (planet)
    case ('mars')
      call fiona_output('h0', 'Ls', curr_time%solar_longitude() * deg)
    end select

    if (advection) then
      call history_write_h0_adv(blocks, itime)
    else if (hydrostatic) then
      call history_write_h0_hydrostatic(blocks, itime)
    else if (nonhydrostatic) then
      call history_write_h0_nonhydrostatic(blocks, itime)
    else
      call history_write_h0_swm(blocks, itime)
    end if

    call fiona_end_output('h0', keep_dataset=.true.)

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write state cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_h0

  subroutine history_write_h1_swm(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
    call fiona_output('h1', 'dmf'     ,  dtend%dmf     (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke      (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dgzdt'   ,  dtend%dgz     (is:ie,js:je,1), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    start = [is,js]
    count = [mesh%half_nlon,mesh%full_nlat]
    call fiona_output('h1', 'dudt   ' ,  dtend%du       (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' , dstate%mfx_lon  (is:ie,js:je,1), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    start = [is,js]
    count = [mesh%full_nlon,mesh%half_nlat]
    call fiona_output('h1', 'dvdt'    ,  dtend%dv       (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' , dstate%mfy_lat  (is:ie,js:je,1), start=start, count=count)
    end associate

    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_swm

  subroutine history_write_h1_hydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
               pstate => blocks(1)%pstate       , &
               dtend  => blocks(1)%dtend (itime), &
               aux    => blocks(1)%aux          , &
               static => blocks(1)%static)
    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dmf'     ,  dtend%dmf      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg'     , dstate%dmg      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmgsdt'  ,  dtend%dmgs     (is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  dtend%dpt      (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%half_nlat,mesh%full_nlev]

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dudt   ' ,  dtend%du       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' , dstate%mfx_lon  (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dvdt'    ,  dtend%dv       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' , dstate%mfy_lat  (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%half_kds; ke = mesh%half_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
    call fiona_output('h1', 'we_lev', dstate%we_lev(is:ie,js:je,ks:ke), start=start, count=count)

    call fiona_output('h1', 'n2', pstate%n2_lev, start=start, count=count)
    call fiona_output('h1', 'ri', pstate%ri_lev, start=start, count=count)

    if (physics_suite /= 'none') then
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('h1', 'dudt_phys' , dtend%dudt_phys (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dvdt_phys' , dtend%dvdt_phys (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dtdt_phys' , dtend%dtdt_phys (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dshdt_phys', dtend%dshdt_phys(is:ie,js:je,ks:ke), start=start, count=count)
    end if
    end associate
    
    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_hydrostatic

  subroutine history_write_h1_nonhydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
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
    call fiona_output('h1', 'dmf'     ,  dtend%dmf      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      ,    aux%ke       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmg'     , dstate%dmg      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmgsdt'  ,  dtend%dmgs     (is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  dtend%dpt      (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_ids; ie = mesh%half_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dudt   ' ,  dtend%du       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' , dstate%mfx_lon  (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%half_jds; je = mesh%half_jde
    ks = mesh%full_kds; ke = mesh%full_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
    call fiona_output('h1', 'dvdt'    ,  dtend%dv       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' , dstate%mfy_lat  (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_ids; ie = mesh%full_ide
    js = mesh%full_jds; je = mesh%full_jde
    ks = mesh%half_kds; ke = mesh%half_kde
    start = [is,js,ks]
    count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
    call fiona_output('h1', 'we_lev', dstate%we_lev(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_gz',  dtend%adv_gz(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_w' ,  dtend%adv_w (is:ie,js:je,ks:ke), start=start, count=count)
    end associate

    call fiona_end_output('h1', keep_dataset=.true.)

  end subroutine history_write_h1_nonhydrostatic

  subroutine history_write_h1(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    if (hydrostatic) then
      call history_write_h1_hydrostatic(blocks, itime)
    else if (nonhydrostatic) then
      call history_write_h1_nonhydrostatic(blocks, itime)
    else
      call history_write_h1_swm(blocks, itime)
    end if

  end subroutine history_write_h1

end module history_mod
