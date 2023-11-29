! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module initial_mod

  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod
  use formula_mod
  use time_mod
  use block_mod
  use tracer_mod
  use latlon_parallel_mod
  use operators_mod
  use process_mod, only: proc

  implicit none

  private

  public initial_write
  public initial_read

  character(30) :: name_style = 'cam'

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
      select case (name_style)
      case ('gmcore')
        call fiona_add_var('i0', 'pt'         , long_name='Potential temperature'                     , units='K'       , dtype=output_i0_dtype, dim_names=cell_dims)
        call fiona_add_var('i0', 'qv'         , long_name='Water vapor dry mixing ratio'              , units='kg kg-1' , dtype=output_i0_dtype, dim_names=cell_dims)
        call fiona_add_var('i0', 'mgs'        , long_name='Surface dry-air weight'                    , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
        call fiona_add_var('i0', 'u'          , long_name='U wind component'                          , units='m s-1'   , dtype=output_i0_dtype, dim_names=lon_dims)
        call fiona_add_var('i0', 'v'          , long_name='V wind component'                          , units='m s-1'   , dtype=output_i0_dtype, dim_names=lat_dims)
        call fiona_add_var('i0', 'z'          , long_name='Height'                                    , units='m'       , dtype=output_i0_dtype, dim_names=lev_dims)
        call fiona_add_var('i0', 'ref_ps_smth', long_name='Smoothed reference surface pressure'       , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
        call fiona_add_var('i0', 'ref_ps_perb', long_name='Perturbed reference surface pressure'      , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      case ('cam')
        call fiona_add_var('i0', 'T'          , long_name='Temperature'                               , units='K'       , dtype=output_i0_dtype, dim_names=cell_dims)
        call fiona_add_var('i0', 'Q'          , long_name='Specific humidity'                         , units='kg/kg'   , dtype=output_i0_dtype, dim_names=cell_dims)
        call fiona_add_var('i0', 'PS'         , long_name='Surface pressure'                          , units='Pa'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
        call fiona_add_var('i0', 'US'         , long_name='Zonal wind, staggered'                     , units='m/s'     , dtype=output_i0_dtype, dim_names=lon_dims)
        call fiona_add_var('i0', 'VS'         , long_name='Meridional wind, staggered'                , units='m/s'     , dtype=output_i0_dtype, dim_names=lat_dims)
        if (idx_qc > 0) &
          call fiona_add_var('i0', 'CLDLIQ'   , long_name='Grid box averaged cloud liquid amount'     , units='1/m3'    , dtype=output_i0_dtype, dim_names=cell_dims)
        if (idx_qi > 0) &
          call fiona_add_var('i0', 'CLDICE'   , long_name='Grid box averaged cloud ice amount'        , units='1/m3'    , dtype=output_i0_dtype, dim_names=cell_dims)
        if (idx_nc > 0) &
          call fiona_add_var('i0', 'NUMLIQ'   , long_name='Grid box averaged cloud liquid number'     , units='1/m3'    , dtype=output_i0_dtype, dim_names=cell_dims)
        if (idx_ni > 0) &
          call fiona_add_var('i0', 'NUMICE'   , long_name='Grid box averaged cloud ice number'        , units='1/m3'    , dtype=output_i0_dtype, dim_names=cell_dims)
      end select
    end if
    select case (name_style)
    case ('gmcore')
      call fiona_add_var('i0', 'zs'           , long_name='Surface height'                            , units='m'       , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'zs_std'       , long_name='Surface height subgrid standard deviation' , units='m2'      , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'dzsdx'        , long_name='Zonal zs gradient'                         , units=''        , dtype=output_i0_dtype, dim_names=lon_dims(1:2))
      call fiona_add_var('i0', 'dzsdy'        , long_name='Meridional zs gradient'                    , units=''        , dtype=output_i0_dtype, dim_names=lat_dims(1:2))
      call fiona_add_var('i0', 'landmask'     , long_name='Land mask'                                 , units='1'       , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    case ('cam')
      call fiona_add_var('i0', 'PHIS'         , long_name='Surface geopotential'                      , units='m2/s2'   , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    end select

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
      associate (mesh     => blocks(iblk)%mesh            , &
                 gzs      => blocks(iblk)%static%gzs      , &
                 zs_std   => blocks(iblk)%static%zs_std   , &
                 landmask => blocks(iblk)%static%landmask , &
                 dzsdx    => blocks(iblk)%static%dzsdx    , &
                 dzsdy    => blocks(iblk)%static%dzsdy    , &
                 mgs      => blocks(iblk)%dstate(1)%mgs   , &
                 u_lon    => blocks(iblk)%dstate(1)%u_lon , &
                 v_lat    => blocks(iblk)%dstate(1)%v_lat , &
                 pt       => blocks(iblk)%dstate(1)%pt    , &
                 t        => blocks(iblk)%dstate(1)%t     , &
                 gz_lev   => blocks(iblk)%dstate(1)%gz_lev, &
                 q        => tracers(iblk)%q              , &
                 qm       => tracers(iblk)%qm             )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

      if (baroclinic) then
        select case (name_style)
        case ('gmcore')
          call fiona_output('i0', 'pt', pt(is:ie,js:je,ks:ke), start=start, count=count)
        case ('cam')
          call fiona_output('i0', 'T' , t (is:ie,js:je,ks:ke), start=start, count=count)
        end select
        if (idx_qv > 0) then
          select case (name_style)
          case ('gmcore')
            call fiona_output('i0', 'qv', q(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
          case ('cam')
            call fiona_output('i0', 'Q' , wet_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qv), qm(is:ie,js:je,ks:ke)), start=start, count=count)
          end select
        end if
        if (idx_qc > 0 .and. fiona_has_var('i0', 'qc')) then
          select case (name_style)
          case ('gmcore')
            call fiona_output('i0', 'qc'    , q(is:ie,js:je,ks:ke,idx_qc), start=start, count=count)
          case ('cam')
            call fiona_output('i0', 'CLDLIQ', wet_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qc), qm(is:ie,js:je,ks:ke)), start=start, count=count)
          end select
        end if
        if (idx_qi > 0) then
          select case (name_style)
          case ('gmcore')
            call fiona_output('i0', 'qi'    , q(is:ie,js:je,ks:ke,idx_qi), start=start, count=count)
          case ('cam')
            call fiona_output('i0', 'CLDICE', wet_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qi), qm(is:ie,js:je,ks:ke)), start=start, count=count)
          end select
        end if
        if (idx_nc > 0) then
          select case (name_style)
          case ('gmcore')
            call fiona_output('i0', 'nc'    , q(is:ie,js:je,ks:ke,idx_nc), start=start, count=count)
          case ('cam')
            call fiona_output('i0', 'NUMLIQ', q(is:ie,js:je,ks:ke,idx_nc), start=start, count=count)
          end select
        end if
        if (idx_ni > 0) then
          select case (name_style)
          case ('gmcore')
            call fiona_output('i0', 'ni'    , q(is:ie,js:je,ks:ke,idx_ni), start=start, count=count)
          case ('cam')
            call fiona_output('i0', 'NUMICE', q(is:ie,js:je,ks:ke,idx_ni), start=start, count=count)
          end select
        end if
        select case (name_style)
        case ('gmcore')
          call fiona_output('i0', 'mgs', mgs(is:ie,js:je), start=start, count=count)
        case ('cam')
          call fiona_output('i0', 'PS' , mgs(is:ie,js:je), start=start, count=count)
        end select
      end if
      select case (name_style)
      case ('gmcore')
        call fiona_output('i0', 'zs'      , gzs     (is:ie,js:je) / g, start=start, count=count)
        call fiona_output('i0', 'zs_std'  , zs_std  (is:ie,js:je)    , start=start, count=count)
        call fiona_output('i0', 'landmask', landmask(is:ie,js:je)    , start=start, count=count)
      case ('cam')
        call fiona_output('i0', 'PHIS'    , gzs     (is:ie,js:je)    , start=start, count=count)
      end select

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_output('i0', 'dzsdx', dzsdx(is:ie,js:je)      , start=start, count=count)
        call fiona_output('i0', 'u_lon', u_lon(is:ie,js:je,ks:ke), start=start, count=count)
      case ('cam')
        call fiona_output('i0', 'US'   , u_lon(is:ie,js:je,ks:ke), start=start, count=count)
      end select

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_output('i0', 'dzsdy', dzsdy(is:ie,js:je)      , start=start, count=count)
        call fiona_output('i0', 'v_lat', v_lat(is:ie,js:je,ks:ke), start=start, count=count)
      case ('cam')
        call fiona_output('i0', 'VS'   , v_lat(is:ie,js:je,ks:ke), start=start, count=count)
      end select

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_output('i0', 'z', gz_lev(is:ie,js:je,ks:ke) / g, start=start, count=count)
      end select
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

    integer iblk, is, ie, js, je, ks, ke, i, j, k
    integer start(3), count(3)
    real(8) time1, time2

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
      associate (block    => blocks(iblk)                 , &
                 mesh     => blocks(iblk)%mesh            , &
                 gzs      => blocks(iblk)%static%gzs      , &
                 zs_std   => blocks(iblk)%static%zs_std   , &
                 landmask => blocks(iblk)%static%landmask , &
                 dzsdx    => blocks(iblk)%static%dzsdx    , &
                 dzsdy    => blocks(iblk)%static%dzsdy    , &
                 mgs      => blocks(iblk)%dstate(1)%mgs   , &
                 ph       => blocks(iblk)%dstate(1)%ph    , &
                 u_lon    => blocks(iblk)%dstate(1)%u_lon , &
                 v_lat    => blocks(iblk)%dstate(1)%v_lat , &
                 pt       => blocks(iblk)%dstate(1)%pt    , &
                 t        => blocks(iblk)%dstate(1)%t     , &
                 gz       => blocks(iblk)%dstate(1)%gz    , &
                 gz_lev   => blocks(iblk)%dstate(1)%gz_lev, &
                 q        => tracers(iblk)%q              , &
                 qm       => tracers(iblk)%qm             )
      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%full_nlat,mesh%full_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_input('i0', 'zs', gzs(is:ie,js:je), start=start, count=count); gzs = gzs * g
        call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)
        call fiona_input('i0', 'zs_std', zs_std(is:ie,js:je), start=start, count=count)
        call fill_halo(block%halo, zs_std, full_lon=.true., full_lat=.true.)
        call fiona_input('i0', 'landmask', landmask(is:ie,js:je), start=start, count=count)
        call fill_halo(block%halo, landmask, full_lon=.true., full_lat=.true.)
      case ('cam')
        call fiona_input('i0', 'PHIS', gzs(is:ie,js:je), start=start, count=count)
        call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)
      end select
      if (baroclinic) then
        select case (name_style)
        case ('gmcore')
          call fiona_input('i0', 'mgs', mgs(is:ie,js:je      ), start=start, count=count)
          call fill_halo(block%halo, mgs, full_lon=.true., full_lat=.true.)
          call fiona_input('i0', 'pt' , pt (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          if (idx_qv > 0) then
            call fiona_input('i0', 'qv' , q(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
            call fill_halo(block%filter_halo, q(:,:,:,idx_qv), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
        case ('cam')
          call fiona_input('i0', 'PS', mgs(is:ie,js:je      ), start=start, count=count)
          call fill_halo(block%halo, mgs, full_lon=.true., full_lat=.true.)
          call fiona_input('i0', 'T' , t  (is:ie,js:je,ks:ke), start=start, count=count)
          if (idx_qv > 0) then
            call fiona_input('i0', 'Q' , q(is:ie,js:je,ks:ke,idx_qv), start=start, count=count)
          end if
          if (idx_qc > 0 .and. fiona_has_var('i0', 'CLDLIQ')) then
            call fiona_input('i0', 'CLDLIQ', q(is:ie,js:je,ks:ke,idx_qc), start=start, count=count)
          end if
          if (idx_qi > 0 .and. fiona_has_var('i0', 'CLDICE')) then
            call fiona_input('i0', 'CLDICE', q(is:ie,js:je,ks:ke,idx_qi), start=start, count=count)
          end if
          if (idx_nc > 0 .and. fiona_has_var('i0', 'NUMLIQ')) then
            call fiona_input('i0', 'NUMLIQ', q(is:ie,js:je,ks:ke,idx_nc), start=start, count=count)
            call fill_halo(block%filter_halo, q(:,:,:,idx_nc), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
          if (idx_ni > 0 .and. fiona_has_var('i0', 'NUMICE')) then
            call fiona_input('i0', 'NUMICE', q(is:ie,js:je,ks:ke,idx_ni), start=start, count=count)
            call fill_halo(block%filter_halo, q(:,:,:,idx_ni), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
          call tracer_calc_qm(block)
          if (idx_qv > 0) then
            q(is:ie,js:je,ks:ke,idx_qv) = dry_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qv), qm(is:ie,js:je,ks:ke))
            call fill_halo(block%filter_halo, q(:,:,:,idx_qv), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
          if (idx_qc > 0) then
            q(is:ie,js:je,ks:ke,idx_qc) = dry_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qc), qm(is:ie,js:je,ks:ke))
            call fill_halo(block%filter_halo, q(:,:,:,idx_qc), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
          if (idx_qi > 0) then
            q(is:ie,js:je,ks:ke,idx_qi) = dry_mixing_ratio(q(is:ie,js:je,ks:ke,idx_qi), qm(is:ie,js:je,ks:ke))
            call fill_halo(block%filter_halo, q(:,:,:,idx_qi), full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
          end if
          call calc_mg (block, block%dstate(1))
          call calc_dmg(block, block%dstate(1))
          call calc_ph (block, block%dstate(1))
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                pt(i,j,k) = modified_potential_temperature(t(i,j,k), ph(i,j,k), q(i,j,k,idx_qv))
              end do
            end do
          end do
          call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
        end select
      else
        call fiona_input('i0', 'z' , gz(is:ie,js:je,ks:ke), start=start, count=count); gz = gz * g
        call fill_halo(block%halo, gz, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if

      is = mesh%half_ids; ie = mesh%half_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%full_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_input('i0', 'u_lon', u_lon(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      case ('cam')
        call fiona_input('i0', 'US'   , u_lon(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      end select

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%half_jds; je = mesh%half_jde
      ks = mesh%full_kds; ke = mesh%full_kde
      start = [is,js,ks]
      count = [mesh%full_nlon,mesh%half_nlat,mesh%full_nlev]

      select case (name_style)
      case ('gmcore')
        call fiona_input('i0', 'v_lat', v_lat(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
      case ('cam')
        call fiona_input('i0', 'VS'   , v_lat(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block%halo, v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
      end select

      is = mesh%full_ids; ie = mesh%full_ide
      js = mesh%full_jds; je = mesh%full_jde
      ks = mesh%half_kds; ke = mesh%half_kde
      start = [is,js,ks]
      count = [mesh%half_nlon,mesh%full_nlat,mesh%half_nlev]

      if (nonhydrostatic) then
        call fiona_input('i0', 'z', gz_lev(is:ie,js:je,ks:ke), start=start, count=count); gz_lev = gz_lev * g
        call fill_halo(block%halo, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
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
