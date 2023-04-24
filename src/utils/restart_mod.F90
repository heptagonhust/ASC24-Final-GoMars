module restart_mod

  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use tracer_mod
  use parallel_mod

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
    real(8) time1, time2
    real(r8), pointer, dimension(:,:,:) :: qv

    if (proc%is_root()) then
      call log_notice('Write restart.')
      call MPI_WTIME(time1, ierr)
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
  if (baroclinic) then
    call fiona_add_dim('r0', 'lev' , size=global_mesh%full_nlev, add_var=.true.)
    call fiona_add_dim('r0', 'ilev', size=global_mesh%half_nlev, add_var=.true.)
    call fiona_add_var('r0', 'u'   , long_name='u wind component'            , units='m s-1' , dim_names=lon_dims_3d , dtype='r8')
    call fiona_add_var('r0', 'v'   , long_name='v wind component'            , units='m s-1' , dim_names=lat_dims_3d , dtype='r8')
    call fiona_add_var('r0', 'mgs' , long_name='surface dry-air weight'      , units='Pa'    , dim_names=cell_dims_2d, dtype='r8')
    call fiona_add_var('r0', 'pt'  , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d, dtype='r8')
  if (nonhydrostatic) then
    call fiona_add_var('r0', 'gz_lev', long_name='geopotential height'       , units='m2 s-2', dim_names=lev_dims_3d , dtype='r8')
    call fiona_add_var('r0', 'w'     , long_name='vertical velocity'         , units='m s-1' , dim_names=lev_dims_3d , dtype='r8')
  end if
  else
    call fiona_add_var('r0', 'u'   , long_name='u wind component'            , units='m s-1' , dim_names=lon_dims_2d , dtype='r8')
    call fiona_add_var('r0', 'v'   , long_name='v wind component'            , units='m s-1' , dim_names=lat_dims_2d , dtype='r8')
    call fiona_add_var('r0', 'gz'  , long_name='geopotential height'         , units='m2 s-2', dim_names=cell_dims_2d, dtype='r8')
  end if
    call fiona_add_var('r0', 'gzs' , long_name='surface geopotential height' , units='m2 s-2', dim_names=cell_dims_2d, dtype='r8')

  if (idx_qv > 0) then
    call fiona_add_var('r0', 'qv'  , long_name='water vapor mixing ratio'    , units='kg kg-1', dim_names=cell_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'moist_u0' , long_name='', units='', dim_names= lon_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'moist_v0' , long_name='', units='', dim_names= lat_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'moist_we0', long_name='', units='', dim_names= lev_dims_3d, dtype='r8')
    call fiona_add_var('r0', 'moist_m0' , long_name='', units='', dim_names=cell_dims_3d, dtype='r8')
  end if

    call fiona_start_output('r0', dble(elapsed_seconds), new_file=.true.)
    call fiona_output('r0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('r0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('r0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('r0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh         , &
                 dstate => blocks(iblk)%dstate(itime), &
                 static => blocks(iblk)%static)

        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

        call fiona_output('r0', 'gzs'   , static%gzs   (is:ie,js:je      ), start=start, count=count)
      if (baroclinic) then
        call fiona_output('r0', 'mgs'   , dstate%mgs   (is:ie,js:je      ), start=start, count=count)
        call fiona_output('r0', 'pt'    , dstate%pt    (is:ie,js:je,ks:ke), start=start, count=count)
      else
        call fiona_output('r0', 'gz'    , dstate%gz    (is:ie,js:je,ks:ke), start=start, count=count)
      end if
      if (idx_qv > 0) then
        call tracer_get_array(iblk, idx_qv, qv)
        call fiona_output('r0', 'qv'    , qv           (is:ie,js:je,ks:ke), start=start, count=count)
        associate (adv_batch => blocks(iblk)%adv_batches(1))
        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
        call fiona_output('r0', 'moist_u0', adv_batch%u0(is:ie,js:je,ks:ke), start=start, count=count)
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
        call fiona_output('r0', 'moist_v0', adv_batch%v0(is:ie,js:je,ks:ke), start=start, count=count)
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
        call fiona_output('r0', 'moist_we0', adv_batch%we0(is:ie,js:je,ks:ke), start=start, count=count)
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
        call fiona_output('r0', 'moist_m0', adv_batch%m0(is:ie,js:je,ks:ke), start=start, count=count)
        end associate
      end if

        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

        call fiona_output('r0', 'u'     , dstate%u_lon(is:ie,js:je,ks:ke), start=start, count=count)

        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

        call fiona_output('r0', 'v'     , dstate%v_lat(is:ie,js:je,ks:ke), start=start, count=count)

        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]

      if (nonhydrostatic) then
        call fiona_output('r0', 'gz_lev', dstate%gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('r0', 'w_lev' , dstate%w_lev (is:ie,js:je,ks:ke), start=start, count=count)
      end if
      end associate
    end do
    call fiona_end_output('r0')
    call process_barrier()
    if (proc%is_root()) then
      call MPI_WTIME(time2, ierr)
      call log_notice('Done write restart cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_write

  subroutine restart_read()

    type(block_type), pointer :: block
    type(mesh_type), pointer :: mesh
    type(dstate_type), pointer :: dstate
    type(static_type), pointer :: static
    type(datetime_type) time
    integer iblk, time_step, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(r8) time_value, time1, time2
    real(r8), pointer, dimension(:,:,:) :: qv
    character(50) time_units

    if (restart_file == 'N/A') then
      call log_error('Parameter restart_file is needed to restart!')
    end if

    if (proc%is_root()) then
      call log_notice('Read restart file ' // trim(restart_file) // '.')
      call cpu_time(time1)
    end if

    call fiona_open_dataset('r0', file_path=restart_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_start_input('r0')

    time_step = 1

    time_value = 0
    call fiona_input('r0', 'time', time_value, time_step=time_step)
    call fiona_get_att('r0', 'time', 'units', time_units)
    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)                     , &
                 mesh   => blocks(iblk)%mesh                , &
                 dstate => blocks(iblk)%dstate(old_time_idx), &
                 static => blocks(iblk)%static)
        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

        call fiona_input('r0', 'gzs', static%gzs(is:ie,js:je), start=start, count=count, time_step=time_step)
        call fill_halo(block%filter_halo, static%gzs, full_lon=.true., full_lat=.true.)
        if (baroclinic) then
          call fiona_input('r0', 'mgs', dstate%mgs(is:ie,js:je      ), start=start, count=count, time_step=time_step)
          call fill_halo(block%halo, dstate%mgs, full_lon=.true., full_lat=.true.)
          call fiona_input('r0', 'pt' , dstate%pt (is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%filter_halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
        else
          call fiona_input('r0', 'gz' , dstate%gz (is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%halo, dstate%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
        end if

        if (idx_qv > 0) then
          call tracer_get_array(iblk, idx_qv, qv)
          call fiona_input('r0', 'qv' , qv(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true.)
          associate (adv_batch => block%adv_batches(1))
          is = mesh%half_ids; ie = mesh%half_ide
          js = mesh%full_jds; je = mesh%full_jde
          ks = mesh%full_kds; ke = mesh%full_kde
          start = [is,js,ks]
          count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]
          call fiona_input('r0', 'moist_u0', adv_batch%u0(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%halo, adv_batch%u0, full_lon=.false., full_lat=.true., full_lev=.true.)
          is = mesh%full_ids; ie = mesh%full_ide
          js = mesh%half_jds; je = mesh%half_jde
          ks = mesh%full_kds; ke = mesh%full_kde
          start = [is,js,ks]
          count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]
          call fiona_input('r0', 'moist_v0', adv_batch%v0(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%halo, adv_batch%v0, full_lon=.true., full_lat=.false., full_lev=.true.)
          is = mesh%full_ids; ie = mesh%full_ide
          js = mesh%full_jds; je = mesh%full_jde
          ks = mesh%half_kds; ke = mesh%half_kde
          start = [is,js,ks]
          count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]
          call fiona_input('r0', 'moist_we0', adv_batch%we0(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
          call fill_halo(block%halo, adv_batch%we0, full_lon=.true., full_lat=.true., full_lev=.false.)
          is = mesh%full_ids; ie = mesh%full_ide
          js = mesh%full_jds; je = mesh%full_jde
          ks = mesh%full_kds; ke = mesh%full_kde
          start = [is,js,ks]
          count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]
          call fiona_input('r0', 'moist_m0', adv_batch%m0(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%halo, adv_batch%m0, full_lon=.true., full_lat=.true., full_lev=.true.)
          adv_batch%uv_step = -1
          adv_batch%mf_step = -1
          adv_batch%we_step = -1
          end associate
        end if

        is = mesh%half_ids; ie = mesh%half_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

        call fiona_input('r0', 'u'  , dstate%u_lon(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block%halo, dstate%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%half_jds; je = mesh%half_jde
        ks = mesh%full_kds; ke = mesh%full_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

        call fiona_input('r0', 'v'  , dstate%v_lat(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block%halo, dstate%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

        is = mesh%full_ids; ie = mesh%full_ide
        js = mesh%full_jds; je = mesh%full_jde
        ks = mesh%half_kds; ke = mesh%half_kde
        start = [is,js,ks]
        count = [mesh%full_nlon,mesh%full_nlat,mesh%half_nlev]

      if (nonhydrostatic) then
        call fiona_input('r0', 'gz_lev', dstate%gz_lev(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block%halo, dstate%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

        call fiona_input('r0', 'w_lev' , dstate%w_lev (is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block%halo, dstate%w_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if
      end associate
    end do
    call fiona_end_input('r0')

    call time_fast_forward(time_value, time_units)
    call process_barrier()
    if (proc%is_root()) then
      call cpu_time(time2)
      call log_notice('Restart to ' // trim(curr_time_str) // ' cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine restart_read

end module restart_mod
