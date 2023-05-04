module topo_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use parallel_mod
  use filter_mod
  use laplace_damp_mod

  implicit none

  private

  public topo_read
  public topo_regrid
  public topo_smooth
  public topo_final

  real(r8), allocatable :: topo_lon(:)   ! Longitude (degree)
  real(r8), allocatable :: topo_lat(:)   ! Latitude (degree)
  real(r8), allocatable :: topo_gzs(:,:) ! Geopotential (m2 s-2)

contains

  subroutine topo_read(min_lon, max_lon, min_lat, max_lat)

    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat

    call topo_final()

    if (proc%is_root()) call log_notice('Use ' // trim(topo_file) // ' as topography.')

    call fiona_open_dataset('topo', file_path=topo_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    select case (topo_type)
    case ('etopo1')
      if (planet /= 'earth') call log_error('Topography file ' // trim(topo_file) // ' is used for the Earth!')
      call fiona_set_dim('topo', 'x', span=[-180, 180], cyclic=.true.)
      call fiona_set_dim('topo', 'y', span=[-90, 90])
      call fiona_start_input('topo')
      call fiona_input_range('topo', 'x', topo_lon, coord_range=[min_lon,max_lon])
      call fiona_input_range('topo', 'y', topo_lat, coord_range=[min_lat,max_lat])
      call fiona_input_range('topo', 'z', topo_gzs, coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
      topo_gzs = topo_gzs * g
    case ('gmted')
      if (planet /= 'earth') call log_error('Topography file ' // trim(topo_file) // ' is used for the Earth!')
      call fiona_set_dim('topo', 'lon', span=[-180, 180], cyclic=.true.)
      call fiona_set_dim('topo', 'lat', span=[-90, 90])
      call fiona_start_input('topo')
      call fiona_input_range('topo', 'lon', topo_lon, coord_range=[min_lon,max_lon])
      call fiona_input_range('topo', 'lat', topo_lat, coord_range=[min_lat,max_lat])
      call fiona_input_range('topo', 'htopo', topo_gzs, coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
      topo_gzs = topo_gzs * g
    case ('mola32')
      if (planet /= 'mars') call log_error('Topography file ' // trim(topo_file) // ' is used for the Mars!')
      call fiona_set_dim('topo', 'longitude', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('topo', 'latitude', span=[-90, 90])
      call fiona_start_input('topo')
      call fiona_input_range('topo', 'longitude', topo_lon, coord_range=[min_lon,max_lon])
      call fiona_input_range('topo', 'latitude', topo_lat, coord_range=[min_lat,max_lat])
      call fiona_input_range('topo', 'alt', topo_gzs, coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
      topo_gzs = topo_gzs * g
    case default
      call log_error('Unknown topo_type "' // trim(topo_type) // '"!', pid=proc%id)
    end select
    call fiona_end_input('topo')

  end subroutine topo_read

  subroutine fill_grid(lon1, lon2, lat1, lat2, gzs, std, lnd, cnt)

    real(r8), intent(in) :: lon1
    real(r8), intent(in) :: lon2
    real(r8), intent(in) :: lat1
    real(r8), intent(in) :: lat2
    real(r8), intent(inout) :: gzs
    real(r8), intent(inout) :: std
    real(r8), intent(inout) :: lnd
    integer , intent(inout) :: cnt

    integer is, ie, js, je, i, j

    do is = 1, size(topo_lon)
      if (lon1 <= topo_lon(is)) exit
    end do
    do ie = size(topo_lon), 1, -1
      if (lon2 >= topo_lon(ie)) exit
    end do
    do js = 1, size(topo_lat)
      if (lat1 <= topo_lat(js)) exit
    end do
    do je = size(topo_lat), 1, -1
      if (lat2 >= topo_lat(je)) exit
    end do

    select case (planet)
    case ('earth')
      do j = js, je
        do i = is, ie
          if (topo_lon(i) < lon1 .or. topo_lon(i) > lon2 .or. topo_lat(j) < lat1 .or. topo_lat(j) > lat2) then
            stop 999
          end if
          if (topo_gzs(i,j) > 0) then
            gzs = gzs + topo_gzs(i,j)
            std = std + topo_gzs(i,j)**2
          end if
        end do
      end do
      lnd = lnd + count(topo_gzs(is:ie,js:je) > 0)
    case ('mars')
      gzs = gzs + sum(topo_gzs(is:ie,js:je))
      std = std + sum(topo_gzs(is:ie,js:je)**2)
      lnd = lnd + (ie - is + 1) * (je - js + 1)
    end select
    cnt = cnt + (ie - is + 1) * (je - js + 1)

  end subroutine fill_grid

  subroutine topo_regrid(block)

    type(block_type), intent(inout) :: block

    real(r8) min_lon, max_lon, min_lat, max_lat
    real(r8) lon1, lon2, lat1, lat2, pole_gzs, pole_std, pole_lnd
    integer i, j, pole_n
    integer n(block%mesh%full_ids:block%mesh%full_ide)

    associate (mesh => block%mesh           , &
               lnd  => block%static%landmask, &
               gzs  => block%static%gzs     , &
               std  => block%static%zs_std)
    lnd = 0; gzs = 0; std = 0
    do j = mesh%full_jds, mesh%full_jde
      lat1 = mesh%half_lat_deg(j-1); lat1 = merge(lat1, -90.0_r8, lat1 /= inf)
      lat2 = mesh%half_lat_deg(j  ); lat2 = merge(lat2,  90.0_r8, lat2 /= inf)
      n = 0
      do i = mesh%full_ids, mesh%full_ide
        lon1 = mesh%half_lon_deg(i-1)
        lon2 = mesh%half_lon_deg(i  )
        call fill_grid(lon1, lon2, lat1, lat2, gzs(i,j), std(i,j), lnd(i,j), n(i))
        if (.not. mesh%is_pole(j)) then
          gzs(i,j) = gzs(i,j) / n(i)
          std(i,j) = (std(i,j) - 2 * gzs(i,j)**2 * n(i) + gzs(i,j)**2) / n(i) / g
          lnd(i,j) = lnd(i,j) / n(i)
        end if
      end do
      if (mesh%is_pole(j)) then
        call zonal_sum(proc%zonal_circle, gzs(mesh%full_ids:mesh%full_ide,j), pole_gzs)
        call zonal_sum(proc%zonal_circle, std(mesh%full_ids:mesh%full_ide,j), pole_std)
        call zonal_sum(proc%zonal_circle, lnd(mesh%full_ids:mesh%full_ide,j), pole_lnd)
        call zonal_sum(proc%zonal_circle, n, pole_n)
        gzs(mesh%full_ids:mesh%full_ide,j) = pole_gzs / pole_n
        std(mesh%full_ids:mesh%full_ide,j) = (pole_std - 2 * pole_gzs**2 * pole_n + pole_gzs**2) / pole_n / g
        lnd(mesh%full_ids:mesh%full_ide,j) = pole_lnd / pole_n
      end if
    end do
    call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)
    call fill_halo(block%halo, std, full_lon=.true., full_lat=.true.)
    call fill_halo(block%halo, lnd, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine topo_regrid

  subroutine topo_smooth(block)

    type(block_type), intent(inout) :: block

    real(r8) wgt
    integer i, j, cyc

    if (proc%is_root()) call log_notice('Filter topography.')

    associate (mesh  => block%mesh        , &
               gzs   => block%static%gzs  , &
               gzs_f => block%dstate(1)%pt)   ! Borrow the array.
    do cyc = 1, topo_smooth_cycles
      call filter_on_cell(block%big_filter, gzs, gzs_f(:,:,1))
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        if (abs(mesh%full_lat_deg(j)) > 60) then
          wgt = sin(pi05 * (1 - (pi05 - abs(mesh%full_lat(j))) / (30 * rad)))
          gzs(:,j) = wgt * gzs_f(:,j,1) + (1 - wgt) * gzs(:,j)
        end if
      end do
      call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)
    end do
    wgt = maxval(gzs / g)
    call global_max(proc%comm, wgt)
    if (proc%is_root()) call log_notice('Maximum zs is ' // to_str(wgt, 10) // '.')
    end associate

  end subroutine topo_smooth

  subroutine topo_final()

    if (allocated(topo_lon)) deallocate(topo_lon)
    if (allocated(topo_lat)) deallocate(topo_lat)
    if (allocated(topo_gzs)) deallocate(topo_gzs)

  end subroutine topo_final

end module topo_mod
