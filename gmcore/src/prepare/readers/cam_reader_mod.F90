module cam_reader_mod

  use fiona
  use flogger
  use const_mod
  use namelist_mod
  use formula_mod
  use process_mod

  implicit none

  integer cam_nlon
  integer cam_nlat
  integer cam_nlev

  real(r8), allocatable, dimension(:    ) :: cam_lon
  real(r8), allocatable, dimension(:    ) :: cam_lat
  real(r8), allocatable, dimension(:    ) :: cam_slon
  real(r8), allocatable, dimension(:    ) :: cam_slat
  real(r8), allocatable, dimension(:    ) :: cam_hyam
  real(r8), allocatable, dimension(:    ) :: cam_hybm
  real(r8), allocatable, dimension(:,:,:) :: cam_us
  real(r8), allocatable, dimension(:,:,:) :: cam_vs
  real(r8), allocatable, dimension(:,:,:) :: cam_t
  real(r8), allocatable, dimension(:,:,:) :: cam_q
  real(r8), allocatable, dimension(:,:,:) :: cam_cldliq
  real(r8), allocatable, dimension(:,:,:) :: cam_cldice
  real(r8), allocatable, dimension(:,:,:) :: cam_numliq
  real(r8), allocatable, dimension(:,:,:) :: cam_numice
  real(r8), allocatable, dimension(:,:  ) :: cam_ps
  real(r8), allocatable, dimension(:,:,:) :: cam_p
  real(r8) cam_p0

contains

  subroutine cam_reader_run(bkg_file, min_lon, max_lon, min_lat, max_lat)

    character(*), intent(in) :: bkg_file
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat

    integer i, j, k, k0
    real(r8) qm1, qm2, qm, tv1, tv2, sum_q_lev, dz, rho
    real(r8), allocatable :: p(:)

    call cam_reader_final()

    if (proc%is_root()) call log_notice('Use CAM initial file ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('cam', file_path=bkg_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('cam', 'lon', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('cam', 'lat', span=[-90, 90])
    call fiona_get_dim('cam', 'lev', size=cam_nlev)
    allocate(cam_hyam(cam_nlev))
    allocate(cam_hybm(cam_nlev))
    call fiona_start_input('cam')
    call fiona_input('cam', 'hyam', cam_hyam)
    call fiona_input('cam', 'hybm', cam_hybm)
    call fiona_input('cam', 'P0'  , cam_p0)
    call fiona_input_range('cam', 'lon' , cam_lon , coord_range=[min_lon, max_lon]); cam_nlon = size(cam_lon)
    call fiona_input_range('cam', 'slon', cam_slon, coord_range=[min_lon, max_lon])
    call fiona_input_range('cam', 'lat' , cam_lat , coord_range=[min_lat, max_lat]); cam_nlat = size(cam_lat)
    call fiona_input_range('cam', 'slat', cam_slat, coord_range=[min_lat, max_lat])
    call fiona_input_range('cam', 'US'  , cam_us  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('cam', 'VS'  , cam_vs  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('cam', 'T '  , cam_t   , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('cam', 'Q '  , cam_q   , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('cam', 'PS'  , cam_ps  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    if (fiona_has_var('cam', 'CLDLIQ')) then
      call fiona_input_range('cam', 'CLDLIQ', cam_cldliq, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    else
      allocate(cam_cldliq(cam_nlon,cam_nlat,cam_nlev)); cam_cldliq = 0
    end if
    if (fiona_has_var('cam', 'CLDICE')) then
      call fiona_input_range('cam', 'CLDICE', cam_cldice, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    else
      allocate(cam_cldice(cam_nlon,cam_nlat,cam_nlev)); cam_cldice = 0
    end if
    if (fiona_has_var('cam', 'NUMLIQ')) then
      call fiona_input_range('cam', 'NUMLIQ', cam_numliq, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    else
      allocate(cam_numliq(cam_nlon,cam_nlat,cam_nlev)); cam_numliq = 0
    end if
    if (fiona_has_var('cam', 'NUMICE')) then
      call fiona_input_range('cam', 'NUMICE', cam_numice, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    else
      allocate(cam_numice(cam_nlon,cam_nlat,cam_nlev)); cam_numice = 0
    end if
    call fiona_end_input('cam')

    ! Change moist mixing ratio to dry mixing ratio.
    do k = 1, cam_nlev
      do j = 1, cam_nlat
        do i = 1, cam_nlon
          qm = cam_q(i,j,k) + cam_cldliq(i,j,k) + cam_cldice(i,j,k)
          cam_q(i,j,k) = cam_q(i,j,k) / (1 - qm)
          cam_cldliq(i,j,k) = cam_cldliq(i,j,k) / (1 - qm)
          cam_cldice(i,j,k) = cam_cldice(i,j,k) / (1 - qm)
        end do
      end do
    end do

    if (proc%is_root()) call log_notice('Calculate CAM dry-air pressure or weight.')
    allocate(cam_p(cam_nlon,cam_nlat,cam_nlev))
    do k = 1, cam_nlev
      cam_p(:,:,k) = cam_hyam(k) * cam_p0 + cam_hybm(k) * cam_ps
    end do

  end subroutine cam_reader_run

  subroutine cam_reader_final()

    if (allocated(cam_lon   )) deallocate(cam_lon   )
    if (allocated(cam_lat   )) deallocate(cam_lat   )
    if (allocated(cam_slon  )) deallocate(cam_slon  )
    if (allocated(cam_slat  )) deallocate(cam_slat  )
    if (allocated(cam_hyam  )) deallocate(cam_hyam  )
    if (allocated(cam_hybm  )) deallocate(cam_hybm  )
    if (allocated(cam_us    )) deallocate(cam_us    )
    if (allocated(cam_vs    )) deallocate(cam_vs    )
    if (allocated(cam_t     )) deallocate(cam_t     )
    if (allocated(cam_q     )) deallocate(cam_q     )
    if (allocated(cam_cldliq)) deallocate(cam_cldliq)
    if (allocated(cam_cldice)) deallocate(cam_cldice)
    if (allocated(cam_numliq)) deallocate(cam_numliq)
    if (allocated(cam_numice)) deallocate(cam_numice)
    if (allocated(cam_ps    )) deallocate(cam_ps    )
    if (allocated(cam_p     )) deallocate(cam_p     )

  end subroutine cam_reader_final

end module cam_reader_mod