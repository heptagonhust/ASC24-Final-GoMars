module mars_data_nasa_mod

  use fiona
  use const_mod
  use namelist_mod
  use parallel_mod

  implicit none

  integer, parameter :: nlev_soil = 40 ! Number of soil layers

  type mars_data_type
    real(r8), allocatable, dimension(:  ) :: alb        ! Surface albedo
    real(r8), allocatable, dimension(:,:) :: tin        ! Thermal interial
    real(r8), allocatable, dimension(:  ) :: z0         ! Surface roughness length (m)
    real(r8), allocatable, dimension(:,:) :: rho_soil   ! Soil density (kg m-3)
    real(r8), allocatable, dimension(:,:) :: cp_soil    ! Soil specific capacity (J kg-1 K-1)
    real(r8), allocatable, dimension(:,:) :: cn_soil    ! Soil thermal conductivity (W m-1 K-1)
  end type mars_data_type

  type(mars_data_type), allocatable :: mars_data(:)

contains

  subroutine mars_data_nasa_init(nblk, ncol, nlev)

    integer, intent(in) :: nblk
    integer, intent(in) :: ncol(:)
    integer, intent(in) :: nlev

    integer iblk

    call mars_data_nasa_final()

    allocate(mars_data(nblk))

    do iblk = 1, nblk
      allocate(mars_data(iblk)%alb     (ncol(iblk)))
      allocate(mars_data(iblk)%tin     (ncol(iblk),nlev_soil))
      allocate(mars_data(iblk)%z0      (ncol(iblk)))
      allocate(mars_data(iblk)%rho_soil(ncol(iblk),nlev_soil))
      allocate(mars_data(iblk)%cp_soil (ncol(iblk),nlev_soil))
      allocate(mars_data(iblk)%cn_soil (ncol(iblk),nlev_soil))
    end do

  end subroutine mars_data_nasa_init

  subroutine mars_data_nasa_final()

    integer iblk

    if (allocated(mars_data)) then
      do iblk = 1, size(mars_data)
        if (allocated(mars_data(iblk)%alb     )) deallocate(mars_data(iblk)%alb     )
        if (allocated(mars_data(iblk)%tin     )) deallocate(mars_data(iblk)%tin     )
        if (allocated(mars_data(iblk)%z0      )) deallocate(mars_data(iblk)%z0      )
        if (allocated(mars_data(iblk)%rho_soil)) deallocate(mars_data(iblk)%rho_soil)
        if (allocated(mars_data(iblk)%cp_soil )) deallocate(mars_data(iblk)%cp_soil )
        if (allocated(mars_data(iblk)%cn_soil )) deallocate(mars_data(iblk)%cn_soil )
      end do
      deallocate(mars_data)
    end if

  end subroutine mars_data_nasa_final

  subroutine mars_data_read()

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

  end subroutine mars_data_read

end module mars_data_nasa_mod