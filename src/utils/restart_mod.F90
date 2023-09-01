! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module writes and reads restart files in NetCDF format.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module restart_mod

  use mpi
  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use tracer_mod
  use latlon_parallel_mod
  use process_mod, only: proc

  implicit none

  private

  public restart_init
  public restart_write
  public restart_read

contains

  subroutine restart_init()

    character(10) time_value, time_units
    real(r8) seconds

    if (restart_interval == 'N/A') then
      if (proc%is_root()) call log_warning('Parameter restart_interval is not set, so no restart file outputted.')
      return
    end if
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(restart_interval, ' ', 1)
    time_units = split_string(restart_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid restart interval ' // trim(restart_interval) // '!')
    end select

    call time_add_alert('restart_write', seconds=seconds)

  end subroutine restart_init

  subroutine restart_write(itime)

    integer, intent(in) :: itime

    integer iblk, is, ie, js, je, ks, ke, ierr
    integer start(3), count(3)
    character(4) lon_dims_3d(4), lat_dims_3d(4), lev_dims_3d(4), cell_dims_3d(4)
    character(4) lon_dims_2d(3), lat_dims_2d(3),                 cell_dims_2d(3)
    character(30) tag
    real(8) time1, time2
    real(r8), pointer, dimension(:,:,:) :: qv

    if (proc%is_root()) then
      call log_notice('Write restart.')
      time1 = MPI_WTIME()
    end if

     lon_dims_3d(1) = 'ilon';  lon_dims_3d(2) =  'lat';  lon_dims_3d(3) =  'lev';  lon_dims_3d(4) = 'time'
     lat_dims_3d(1) =  'lon';  lat_dims_3d(2) = 'ilat';  lat_dims_3d(3) =  'lev';  lat_dims_3d(4) = 'time'
     lev_dims_3d(1) =  'lon';  lev_dims_3d(2) =  'lat';  lev_dims_3d(3) = 'ilev';  lev_dims_3d(4) = 'time'
    cell_dims_3d(1) =  'lon'; cell_dims_3d(2) =  'lat'; cell_dims_3d(3) =  'lev'; cell_dims_3d(4) = 'time'
     lon_dims_2d(1) = 'ilon';  lon_dims_2d(2) =  'lat';  lon_dims_2d(3) = 'time'
     lat_dims_2d(1) =  'lon';  lat_dims_2d(2) = 'ilat';  lat_dims_2d(3) = 'time'
    cell_dims_2d(1) =  'lon'; cell_dims_2d(2) =  'lat'; cell_dims_2d(3) = 'time'

    call fiona_create_dataset('r0', desc=case_desc, file_prefix=trim(case_name) // '.' // trim(curr_time_str), &
      mpi_comm=proc%comm, ngroup=output_ngroup)

    call fiona_add_att('r0', 'time_step_size', dt_dyn)
    call fiona_add_att('r0', 'restart_interval', restart_interval)
    call fiona_add_dim('r0', 'time', add_var=.true.)
    call fiona_add_dim('r0', 'lon' , size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'lat' , size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilon', size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilat', size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    call fiona_add_var('r0', 'time_step', long_name='', units='', dim_names=['time'], dtype='i4')
    if (baroclinic) then
      call fiona_add_dim('r0', 'lev' , size=global_mesh%full_nlev, add_var=.true.)
      call fiona_add_dim('r0', 'ilev', size=global_mesh%half_nlev, add_var=.true.)
      call fiona_add_var('r0', 'u'   , long_name='u wind component'            , units='m s-1' , dim_names=lon_dims_3d , dtype='r8')
      call fiona_add_var('r0', 'v'   , long_name='v wind component'            , units='m s-1' , dim_names=lat_dims_3d , dtype='r8')
      call fiona_add_var('r0', 'mgs' , long_name='surface dry-air weight'      , units='Pa'    , dim_names=cell_dims_2d, dtype='r8')
      call fiona_add_var('r0', 'pt'  , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d, dtype='r8')
      call fiona_add_var('r0', 'pt_old_m', long_name='', units='', dim_names=cell_dims_3d, dtype='r8')
      if (nonhydrostatic) then
        call fiona_add_var('r0', 'gz_lev', long_name='geopotential height'       , units='m2 s-2', dim_names=lev_dims_3d , dtype='r8')
        call fiona_add_var('r0', 'w'     , long_name='vertical velocity'         , units='m s-1' , dim_names=lev_dims_3d , dtype='r8')
      end if
    else
      call fiona_add_var('r0', 'u'   , long_name='u wind component'            , units='m s-1' , dim_names=lon_dims_2d , dtype='r8')
      call fiona_add_var('r0', 'v'   , long_name='v wind component'            , units='m s-1' , dim_names=lat_dims_2d , dtype='r8')
      call fiona_add_var('r0', 'gz'  , long_name='geopotential height'         , units='m2 s-2', dim_names=cell_dims_2d, dtype='r8')
    end if
    call fiona_add_var('r0', 'mfx_lon', long_name='', units='', dim_names= lon_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'mfy_lat', long_name='', units='', dim_names= lat_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'we_lev' , long_name='', units='', dim_names= lev_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'gzs'    , long_name='surface geopotential height' , units='m2 s-2', dim_names=cell_dims_2d, dtype='r8')

    ! FIXME: Support other tracers.
    if (idx_qv > 0) then
      tag = blocks(1)%adv_batches(1)%name
      call fiona_add_var('r0', 'qv'  , long_name='water vapor mixing ratio'    , units='kg kg-1', dim_names=cell_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_u'  , long_name='', units='', dim_names= lon_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_mfx', long_name='', units='', dim_names= lon_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_v'  , long_name='', units='', dim_names= lat_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_mfy', long_name='', units='', dim_names= lat_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_we' , long_name='', units='', dim_names= lev_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_accum_m'  , long_name='', units='', dim_names=cell_dims_3d, dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_uv_step'  , long_name='', units='', dim_names=    ['time'], dtype='i4')
      call fiona_add_var('r0', trim(tag)//'_mf_step'  , long_name='', units='', dim_names=    ['time'], dtype='i4')
      call fiona_add_var('r0', trim(tag)//'_we_step'  , long_name='', units='', dim_names=    ['time'], dtype='i4')
      call fiona_add_var('r0', trim(tag)//'_last_time', long_name='', units='', dim_names=    ['time'], dtype='r8')
      call fiona_add_var('r0', trim(tag)//'_old_m', long_name='', units='', dim_names=cell_dims_3d, dtype='r8')
    end if

    if (physics_suite /= 'N/A') then
      call fiona_add_var('r0', 'phys_last_time' , long_name='', units='', dim_names=    ['time'], dtype='r8')
    end if

    call fiona_start_output('r0', dble(elapsed_seconds), new_file=.true.)
    call fiona_output('r0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('r0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('r0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('r0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    call fiona_output('r0', 'time_step', time_step)
    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 aux    => blocks(iblk)%aux          , &
                 dstate => blocks(iblk)%dstate(itime), &
                 static => blocks(iblk)%static)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('r0', 'gzs'   , static%gzs(is:ie,js:je      ), start=start, count=count)
      if (baroclinic) then
        call fiona_output('r0', 'mgs'   , dstate%mgs(is:ie,js:je      ), start=start, count=count)
        call fiona_output('r0', 'pt'    , dstate%pt (is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('r0', 'pt_old_m', blocks(iblk)%adv_batch_pt%old_m(is:ie,js:je,ks:ke), start=start, count=count)
      else
        call fiona_output('r0', 'gz'    , dstate%gz (is:ie,js:je,ks:ke), start=start, count=count)
      end if
      if (idx_qv > 0) then
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        call tracer_get_array(iblk, idx_qv, qv, __FILE__, __LINE__)
        call fiona_output('r0', 'qv'    , qv        (is:ie,js:je,ks:ke), start=start, count=count)
        associate (adv_batch => blocks(iblk)%adv_batches(1))
        tag = adv_batch%name
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
        if (adv_batch%uv_step == -1) then
          call fiona_output('r0', trim(tag)//'_accum_u'  , adv_batch%u0 (is:ie,js:je,ks:ke), start=start, count=count)
        else if (adv_batch%uv_step >= 1) then
          call fiona_output('r0', trim(tag)//'_accum_u'  , adv_batch%u  (is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output('r0', trim(tag)//'_accum_mfx', adv_batch%mfx(is:ie,js:je,ks:ke), start=start, count=count)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
        if (adv_batch%uv_step == -1) then
          call fiona_output('r0', trim(tag)//'_accum_v'  , adv_batch%v0 (is:ie,js:je,ks:ke), start=start, count=count)
        else if (adv_batch%uv_step >= 1) then
          call fiona_output('r0', trim(tag)//'_accum_v'  , adv_batch%v  (is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output('r0', trim(tag)//'_accum_mfy', adv_batch%mfy(is:ie,js:je,ks:ke), start=start, count=count)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
        if (adv_batch%we_step == -1) then
          call fiona_output('r0', trim(tag)//'_accum_we', adv_batch%we0(is:ie,js:je,ks:ke), start=start, count=count)
        else if (adv_batch%we_step >= 1) then
          call fiona_output('r0', trim(tag)//'_accum_we', adv_batch%we (is:ie,js:je,ks:ke), start=start, count=count)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        if (adv_batch%we_step == -1) then
          call fiona_output('r0', trim(tag)//'_accum_m' , adv_batch%m0 (is:ie,js:je,ks:ke), start=start, count=count)
        else if (adv_batch%we_step >= 1) then
          call fiona_output('r0', trim(tag)//'_accum_m' , adv_batch%m  (is:ie,js:je,ks:ke), start=start, count=count)
        end if
        call fiona_output('r0', trim(tag)//'_old_m', adv_batch%old_m(is:ie,js:je,ks:ke), start=start, count=count)

        call fiona_output('r0', trim(tag)//'_uv_step'  , adv_batch%uv_step)
        call fiona_output('r0', trim(tag)//'_mf_step'  , adv_batch%mf_step)
        call fiona_output('r0', trim(tag)//'_we_step'  , adv_batch%we_step)
        call fiona_output('r0', trim(tag)//'_last_time', time_get_alert_last_time_timestamp(adv_batch%name))
        end associate
      end if

      if (physics_suite /= 'N/A') then
        call fiona_output('r0', 'phys_last_time', time_get_alert_last_time_timestamp('phys'))
      end if

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_output('r0', 'u'      , dstate%u_lon  (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('r0', 'mfx_lon',    aux%mfx_lon(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_output('r0', 'v'      , dstate%v_lat  (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('r0', 'mfy_lat',    aux%mfy_lat(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
      call fiona_output('r0', 'we_lev', dstate%we_lev(is:ie,js:je,ks:ke), start=start, count=count)

      if (nonhydrostatic) then
        call fiona_output('r0', 'gz_lev', dstate%gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('r0', 'w_lev' , dstate%w_lev (is:ie,js:je,ks:ke), start=start, count=count)
      end if
      end associate
    end do
    call fiona_end_output('r0')
    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Done write restart cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_write

  subroutine restart_read()

    type(block_type), pointer :: block
    type(mesh_type), pointer :: mesh
    type(dstate_type), pointer :: dstate
    type(static_type), pointer :: static
    type(datetime_type) time
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(r8) time_value, time1, time2, tmp
    real(r8), pointer, dimension(:,:,:) :: qv
    character(50) time_units, tag

    if (restart_file == 'N/A') then
      call log_error('Parameter restart_file is needed to restart!')
    end if

    if (proc%is_root()) then
      call log_notice('Read restart file ' // trim(restart_file) // '.')
      time1 = MPI_WTIME()
    end if

    call fiona_open_dataset('r0', file_path=restart_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_start_input('r0')

    call fiona_input('r0', 'time', time_value)
    call fiona_get_att('r0', 'time', 'units', time_units)
    call fiona_input('r0', 'time_step', time_step)

    call time_fast_forward(time_value, time_units)

    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)                     , &
                 mesh   => blocks(iblk)%mesh                , &
                 aux    => blocks(iblk)%aux                 , &
                 dstate => blocks(iblk)%dstate(old_time_idx), &
                 static => blocks(iblk)%static)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_input('r0', 'gzs', static%gzs(is:ie,js:je), start=start, count=count)
      call fill_halo(block%filter_halo, static%gzs, full_lon=.true., full_lat=.true.)
      if (baroclinic) then
        call fiona_input('r0', 'mgs', dstate%mgs(is:ie,js:je      ), start=start, count=count)
        call fill_halo(block%halo, dstate%mgs, full_lon=.true., full_lat=.true.)
        call fiona_input('r0', 'pt' , dstate%pt (is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%filter_halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
        call fiona_input('r0', 'pt_old_m', block%adv_batch_pt%old_m(is:ie,js:je,ks:ke), start=start, count=count)
      else
        call fiona_input('r0', 'gz' , dstate%gz (is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, dstate%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if

      if (idx_qv > 0) then
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        call tracer_get_array(iblk, idx_qv, qv, __FILE__, __LINE__)
        call fiona_input('r0', 'qv' , qv(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
        associate (adv_batch => block%adv_batches(1))
        tag = adv_batch%name
        call fiona_input('r0', trim(tag)//'_uv_step', adv_batch%uv_step)
        call fiona_input('r0', trim(tag)//'_mf_step', adv_batch%mf_step)
        call fiona_input('r0', trim(tag)//'_we_step', adv_batch%we_step)
        call fiona_input('r0', trim(tag)//'_last_time', tmp)
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
        if (adv_batch%uv_step == -1) then
          call fiona_input('r0', trim(tag)//'_accum_u'  , adv_batch%u0 (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%u0 , full_lon=.false., full_lat=.true., full_lev=.true.)
        else if (adv_batch%uv_step >= 1) then
          call fiona_input('r0', trim(tag)//'_accum_u'  , adv_batch%u  (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%u  , full_lon=.false., full_lat=.true., full_lev=.true.)
          call fiona_input('r0', trim(tag)//'_accum_mfx', adv_batch%mfx(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%mfx, full_lon=.false., full_lat=.true., full_lev=.true.)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
        if (adv_batch%uv_step == -1) then
          call fiona_input('r0', trim(tag)//'_accum_v'  , adv_batch%v0 (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%v0 , full_lon=.true., full_lat=.false., full_lev=.true.)
        else if (adv_batch%uv_step >= 1) then
          call fiona_input('r0', trim(tag)//'_accum_v'  , adv_batch%v  (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%v  , full_lon=.true., full_lat=.false., full_lev=.true.)
          call fiona_input('r0', trim(tag)//'_accum_mfy', adv_batch%mfy(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%mfy, full_lon=.true., full_lat=.false., full_lev=.true.)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
        if (adv_batch%we_step == -1) then
          call fiona_input('r0', trim(tag)//'_accum_we', adv_batch%we0(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%we0, full_lon=.true., full_lat=.true., full_lev=.false.)
        else if (adv_batch%we_step >= 1) then
          call fiona_input('r0', trim(tag)//'_accum_we', adv_batch%we (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%we , full_lon=.true., full_lat=.true., full_lev=.false.)
        end if
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        if (adv_batch%we_step == -1) then
          call fiona_input('r0', trim(tag)//'_accum_m' , adv_batch%m0 (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%m0 , full_lon=.true., full_lat=.true., full_lev=.true.)
        else if (adv_batch%we_step >= 1) then
          call fiona_input('r0', trim(tag)//'_accum_m' , adv_batch%m  (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%m  , full_lon=.true., full_lat=.true., full_lev=.true.)
        end if
        call fiona_input('r0', trim(tag)//'_old_m', adv_batch%old_m(is:ie,js:je,ks:ke), start=start, count=count)
        call time_set_alert_last_time(adv_batch%name, tmp)
        end associate
      end if

      if (physics_suite /= 'N/A') then
        call fiona_input('r0', 'phys_last_time', tmp)
        call time_set_alert_last_time('phys', tmp)
      end if

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
      call fiona_input('r0', 'u'      , dstate%u_lon  (is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo, dstate%u_lon  , full_lon=.false., full_lat=.true., full_lev=.true.)
      call fiona_input('r0', 'mfx_lon',    aux%mfx_lon(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo,    aux%mfx_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
      call fiona_input('r0', 'v'      , dstate%v_lat  (is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo, dstate%v_lat  , full_lon=.true., full_lat=.false., full_lev=.true.)
      call fiona_input('r0', 'mfy_lat',    aux%mfy_lat(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo,    aux%mfy_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
      call fiona_input('r0', 'we_lev', dstate%we_lev(is:ie,js:je,ks:ke), start=start, count=count)

      if (nonhydrostatic) then
        call fiona_input('r0', 'gz_lev', dstate%gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, dstate%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

        call fiona_input('r0', 'w_lev' , dstate%w_lev (is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, dstate%w_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if
      end associate
    end do
    call fiona_end_input('r0')

    if (proc%is_root()) then
      time2 = MPI_WTIME()
      call log_notice('Restart to ' // trim(curr_time_str) // ' cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_read

end module restart_mod
