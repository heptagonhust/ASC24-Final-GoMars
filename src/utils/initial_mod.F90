module initial_mod

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

  public initial_write
  public initial_read

contains

  subroutine initial_write(initial_file, initial_time)

    character(*), intent(in) :: initial_file
    character(*), intent(in) :: initial_time

    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4)
    character(4) lat_dims(4)
    character(4) lev_dims(4)
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(8) time1, time2
    real(r8), pointer, dimension(:,:,:) :: qv

    cell_dims   (1) =  'lon';  cell_dims   (2) =  'lat';  cell_dims   (3) =  'lev';     cell_dims(4) = 'time'
     lon_dims   (1) = 'ilon';   lon_dims   (2) =  'lat';   lon_dims   (3) =  'lev';      lon_dims(4) = 'time'
     lat_dims   (1) =  'lon';   lat_dims   (2) = 'ilat';   lat_dims   (3) =  'lev';      lat_dims(4) = 'time'
     lev_dims   (1) =  'lon';   lev_dims   (2) =  'lat';   lev_dims   (3) = 'ilev';      lev_dims(4) = 'time'
    cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'

    if (proc%is_root()) then
      call log_notice('Write ' // trim(initial_file) // '.')
      call cpu_time(time1)
    end if

    call fiona_create_dataset('i0', file_path=initial_file, start_time=initial_time, time_units='hours', mpi_comm=proc%comm, ngroup=output_ngroup)
    call fiona_add_dim('i0', 'time', add_var=.true.)
    call fiona_add_dim('i0',  'lon', size=global_mesh%full_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0',  'lat', size=global_mesh%full_nlat, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilon', size=global_mesh%half_nlon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilat', size=global_mesh%half_nlat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('i0',  'lev', size=global_mesh%full_nlev, add_var=.true., decomp=.false.)
      call fiona_add_dim('i0', 'ilev', size=global_mesh%half_nlev, add_var=.true., decomp=.false.)
      call fiona_add_var('i0', 'pt'  , long_name='potential temperature'       , units='K'    , dtype=output_i0_dtype, dim_names=cell_dims)
      call fiona_add_var('i0', 'qv'  , long_name='water vapor mixing ratio'    , units='1'    , dtype=output_i0_dtype, dim_names=cell_dims)
      call fiona_add_var('i0', 'mgs' , long_name='surface dry-air weight'      , units='Pa'   , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'u'   , long_name='u wind component'            , units='m s-1', dtype=output_i0_dtype, dim_names=lon_dims)
      call fiona_add_var('i0', 'v'   , long_name='v wind component'            , units='m s-1', dtype=output_i0_dtype, dim_names=lat_dims)
      call fiona_add_var('i0', 'z'   , long_name='height'                      , units='m'    , dtype=output_i0_dtype, dim_names=lev_dims)
      call fiona_add_var('i0', 'ref_ps_smth', long_name='smoothed reference surface pressure', units='Pa', dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'ref_ps_perb', long_name='perturbed reference surface pressure', units='Pa', dtype=output_i0_dtype, dim_names=cell_dims_2d)
    end if
    call fiona_add_var('i0', 'zs'    , long_name='surface height', units='m' , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'zs_std', long_name='surface height subgrid standard deviation', units='m2' , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'dzsdlon', long_name='zonal zs gradient'     , units='', dim_names=['ilon', 'lat '], dtype=output_i0_dtype)
    call fiona_add_var('i0', 'dzsdlat', long_name='meridional zs gradient', units='', dim_names=['lon ', 'ilat'], dtype=output_i0_dtype)
    call fiona_add_var('i0', 'landmask', long_name='land mask', units='-', dtype=output_i0_dtype, dim_names=cell_dims_2d)

    call fiona_start_output('i0', 0.0d0)
    call fiona_output('i0',  'lon', global_mesh%full_lon_deg(1:global_mesh%full_nlon))
    call fiona_output('i0',  'lat', global_mesh%full_lat_deg(1:global_mesh%full_nlat))
    call fiona_output('i0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%half_nlon))
    call fiona_output('i0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%half_nlat))
    if (baroclinic) then
      call fiona_output('i0',  'lev', global_mesh%full_lev(1:global_mesh%full_nlev))
      call fiona_output('i0', 'ilev', global_mesh%half_lev(1:global_mesh%half_nlev))
    end if

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh     , &
                 dstate => blocks(iblk)%dstate(1), &
                 static => blocks(iblk)%static)
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

      if (baroclinic) then
        call fiona_output('i0', 'pt', dstate%pt (is:ie,js:je,ks:ke), start=start, count=count)
        if (idx_qv > 0) then
          call tracer_get_array(iblk, idx_qv, qv, __FILE__, __LINE__)
          call fiona_output('i0', 'qv', qv(is:ie,js:je,ks:ke), start=start, count=count)
        end if
        call fiona_output('i0', 'mgs', dstate%mgs(is:ie,js:je      ), start=start, count=count)
        call fiona_output('i0', 'ref_ps_smth', static%ref_ps_smth(is:ie,js:je), start=start, count=count)
        call fiona_output('i0', 'ref_ps_perb', static%ref_ps_perb(is:ie,js:je), start=start, count=count)
      end if
      call fiona_output('i0', 'zs'      , static%gzs     (is:ie,js:je) / g, start=start, count=count)
      call fiona_output('i0', 'zs_std'  , static%zs_std  (is:ie,js:je)    , start=start, count=count)
      call fiona_output('i0', 'landmask', static%landmask(is:ie,js:je)    , start=start, count=count)

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      call fiona_output('i0', 'dzsdlon', static%dzsdlon(is:ie,js:je), start=start, count=count)
      call fiona_output('i0', 'u', dstate%u_lon(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      call fiona_output('i0', 'dzsdlat', static%dzsdlat(is:ie,js:je), start=start, count=count)
      call fiona_output('i0', 'v', dstate%v_lat(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      call fiona_output('i0', 'z', dstate%gz_lev(is:ie,js:je,ks:ke) / g, start=start, count=count)
      end associate
    end do
    call fiona_end_output('i0')

    if (proc%is_root()) then
      call cpu_time(time2)
      call log_notice('Done write initial data cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine initial_write

  subroutine initial_read(initial_file_)

    character(*), intent(in), optional :: initial_file_

    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(8) time1, time2
    real(r8), pointer, dimension(:,:,:) :: qv

    if (proc%is_root()) call cpu_time(time1)

    if (present(initial_file_)) then
      call fiona_open_dataset('i0', file_path=initial_file_, mpi_comm=proc%comm, ngroup=input_ngroup)
      if (proc%is_root()) call log_notice('Read initial data from ' // trim(initial_file_) // '.')
    else
      call fiona_open_dataset('i0', file_path=initial_file, mpi_comm=proc%comm, ngroup=input_ngroup)
      if (proc%is_root()) call log_notice('Read initial data from ' // trim(initial_file) // '.')
    end if
    call fiona_start_input('i0')

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

      call fiona_input('i0', 'zs', static%gzs(is:ie,js:je), start=start, count=count)
      static%gzs = static%gzs * g
      call fill_halo(block%filter_halo, static%gzs, full_lon=.true., full_lat=.true.)
      call fiona_input('i0', 'zs_std', static%zs_std(is:ie,js:je), start=start, count=count)
      call fill_halo(block%halo, static%zs_std, full_lon=.true., full_lat=.true.)
      call fiona_input('i0', 'landmask', static%landmask(is:ie,js:je), start=start, count=count)
      call fill_halo(block%halo, static%landmask, full_lon=.true., full_lat=.true.)
      if (baroclinic) then
        call fiona_input('i0', 'mgs', dstate%mgs(is:ie,js:je      ), start=start, count=count)
        call fill_halo(block%halo, dstate%mgs, full_lon=.true., full_lat=.true.)
        call fiona_input('i0', 'pt' , dstate%pt (is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%filter_halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
        if (idx_qv > 0) then
          call tracer_get_array(iblk, idx_qv, qv, __FILE__, __LINE__)
          call fiona_input('i0', 'qv' , qv(is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true.)
        end if
      else
        call fiona_input('i0', 'z' , dstate%gz (is:ie,js:je,ks:ke), start=start, count=count)
        dstate%gz = dstate%gz * g
        call fill_halo(block%halo, dstate%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      call fiona_input('i0', 'u'  , dstate%u_lon(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo, dstate%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      call fiona_input('i0', 'v'  , dstate%v_lat(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block%halo, dstate%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      if (nonhydrostatic) then
        call fiona_input('i0', 'z', dstate%gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
        dstate%gz_lev = dstate%gz_lev * g
        call fill_halo(block%halo, dstate%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if
      end associate
    end do
    call fiona_end_input('i0')

    if (proc%is_root()) then
      call cpu_time(time2)
      call log_notice('Done read initial file cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine initial_read

end module initial_mod
