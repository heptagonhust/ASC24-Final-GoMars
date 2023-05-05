module mars_data_mod

  use fiona
  use const_mod
  use block_mod
  use namelist_mod
  use parallel_mod
  use mars_nasa_mod

  implicit none

  integer, parameter :: nlev_soil = 40 ! Number of soil layers
  real(r8), allocatable, dimension(:) :: soil_dz        ! Soil layer thickness (m)
  real(r8), allocatable, dimension(:) :: soil_z         ! Soil layer depth (m)
  real(r8), allocatable, dimension(:) :: soil_z_lev     ! Soil interface depth (m)

  type mars_data_type
    integer :: ncol = 0 ! Number of columns
    real(r8), pointer    , dimension(:  ) :: lon        ! Longitude (deg)
    real(r8), pointer    , dimension(:  ) :: lat        ! Latitude (deg)
    real(r8), allocatable, dimension(:  ) :: alb        ! Surface albedo
    real(r8), allocatable, dimension(:,:) :: tin        ! Thermal interial
    real(r8), allocatable, dimension(:  ) :: z0         ! Surface roughness length (m)
    real(r8), allocatable, dimension(:,:) :: soil_rho   ! Soil density (kg m-3)
    real(r8), allocatable, dimension(:,:) :: soil_cp    ! Soil specific capacity (J kg-1 K-1)
    real(r8), allocatable, dimension(:,:) :: soil_cn    ! Soil thermal conductivity (W m-1 K-1)
  end type mars_data_type

  type(mars_data_type), allocatable :: mars_data(:)

contains

  subroutine mars_data_init(blocks)

    type(block_type), intent(in), target :: blocks(:)

    integer iblk, ilev

    call mars_data_final()

    allocate(soil_dz   (nlev_soil  ))
    allocate(soil_z    (nlev_soil  ))
    allocate(soil_z_lev(nlev_soil+1))
    soil_dz(1) = factl * skind
    do ilev = 2, nlev_soil
      soil_dz(ilev) = soil_dz(ilev-1) * 1.2_r8
    end do
    soil_z_lev(1) = 0
    do ilev = 2, nlev_soil + 1
      soil_z_lev(ilev) = soil_z_lev(ilev-1) + soil_dz(ilev-1)
    end do
    do ilev = 1, nlev_soil
      soil_z(ilev) = 0.5_r8 * (soil_z_lev(ilev) + soil_z_lev(ilev+1))
    end do

    allocate(mars_data(size(blocks)))

    do iblk = 1, size(blocks)
      associate (ncol => blocks(iblk)%pstate%ncol, &
                 lon  => blocks(iblk)%pstate%lon , &
                 lat  => blocks(iblk)%pstate%lat )
      mars_data(iblk)%ncol = ncol
      mars_data(iblk)%lon => lon
      mars_data(iblk)%lat => lat
      allocate(mars_data(iblk)%alb     (ncol))
      allocate(mars_data(iblk)%tin     (ncol,nlev_soil))
      allocate(mars_data(iblk)%z0      (ncol))
      allocate(mars_data(iblk)%soil_rho(ncol,nlev_soil))
      allocate(mars_data(iblk)%soil_cp (ncol,nlev_soil))
      allocate(mars_data(iblk)%soil_cn (ncol,nlev_soil))
      end associate
    end do

  end subroutine mars_data_init

  subroutine mars_data_final()

    integer iblk

    if (allocated(soil_dz   )) deallocate(soil_dz   )
    if (allocated(soil_z    )) deallocate(soil_z    )
    if (allocated(soil_z_lev)) deallocate(soil_z_lev)
    if (allocated(mars_data)) then
      do iblk = 1, size(mars_data)
        if (allocated(mars_data(iblk)%alb     )) deallocate(mars_data(iblk)%alb     )
        if (allocated(mars_data(iblk)%tin     )) deallocate(mars_data(iblk)%tin     )
        if (allocated(mars_data(iblk)%z0      )) deallocate(mars_data(iblk)%z0      )
        if (allocated(mars_data(iblk)%soil_rho)) deallocate(mars_data(iblk)%soil_rho)
        if (allocated(mars_data(iblk)%soil_cp )) deallocate(mars_data(iblk)%soil_cp )
        if (allocated(mars_data(iblk)%soil_cn )) deallocate(mars_data(iblk)%soil_cn )
      end do
      deallocate(mars_data)
    end if

  end subroutine mars_data_final

  subroutine mars_data_read()

    use latlon_interp_mod

    integer iblk
    integer data_nlon, data_nlat
    real(r8), allocatable :: data_lon(:)
    real(r8), allocatable :: data_lat(:)
    real(r8), allocatable :: data_alb(:,:)
    real(r8), allocatable :: data_tin(:,:)
    real(r8), allocatable :: data_z0 (:,:)

    ! Read in thermal inertial, surface albedo and roughness.
    call fiona_open_dataset('mars', trim(gmcore_data_dir) // 'mars/surface.nc', mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('mars', 'longitude', span=[-180, 180], cyclic=.true.); data_nlon = size(data_lon)
    call fiona_set_dim('mars', 'latitude' , span=[  90, -90],   flip=.true.); data_nlat = size(data_lat)
    call fiona_start_input('mars')
    call fiona_input_range('mars', 'longitude', data_lon, coord_range=[min_lon,max_lon])
    call fiona_input_range('mars', 'latitude' , data_lat, coord_range=[min_lat,max_lat])
    call fiona_input_range('mars', 'albedo'   , data_alb, coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
    call fiona_input_range('mars', 'thermal'  , data_tin, coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
    call fiona_input_range('mars', 'z0'       , data_z0 , coord_range_1=[min_lon,max_lon], coord_range_2=[min_lat,max_lat])
    call fiona_close_dataset('mars')

    do iblk = 1, size(mars_data)
      call latlon_interp_bilinear_column(data_lon, data_lat, data_alb, mars_data(iblk)%lon, mars_data(iblk)%lat, mars_data(iblk)%alb)
      call latlon_interp_bilinear_column(data_lon, data_lat, data_tin, mars_data(iblk)%lon, mars_data(iblk)%lat, mars_data(iblk)%tin(:,1))
      call latlon_interp_bilinear_column(data_lon, data_lat, data_z0 , mars_data(iblk)%lon, mars_data(iblk)%lat, mars_data(iblk)%z0 )
    end do

    deallocate(data_lon, data_lat, data_alb, data_tin, data_z0)

  end subroutine mars_data_read

end module mars_data_mod