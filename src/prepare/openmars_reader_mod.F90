module openmars_reader_mod

  use fiona
  use flogger
  use const_mod
  use block_mod
  use process_mod

  implicit none

  integer openmars_nlon
  integer openmars_nlat
  integer openmars_nlev

  real(r8), allocatable, dimension(:    ) :: openmars_lon
  real(r8), allocatable, dimension(:    ) :: openmars_lat
  real(r8), allocatable, dimension(:    ) :: openmars_lev
  real(r8), allocatable, dimension(:,:,:) :: openmars_u
  real(r8), allocatable, dimension(:,:,:) :: openmars_v
  real(r8), allocatable, dimension(:,:,:) :: openmars_t
  real(r8), allocatable, dimension(:,:,:) :: openmars_p
  real(r8), allocatable, dimension(:,:  ) :: openmars_ps

contains

  subroutine openmars_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer k
    real(r8), allocatable :: tmp(:,:,:)

    call openmars_reader_final()

    if (is_root_proc()) call log_notice('Use OpenMARS ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('openmars', file_path=bkg_file)
    call fiona_get_dim('openmars', 'lon', size=openmars_nlon)
    call fiona_get_dim('openmars', 'lat', size=openmars_nlat)
    call fiona_get_dim('openmars', 'lev', size=openmars_nlev)

    allocate(openmars_lon(openmars_nlon))
    allocate(openmars_lat(openmars_nlat))
    allocate(openmars_lev(openmars_nlev))
    allocate(openmars_u  (openmars_nlon,openmars_nlat,openmars_nlev))
    allocate(openmars_v  (openmars_nlon,openmars_nlat,openmars_nlev))
    allocate(openmars_t  (openmars_nlon,openmars_nlat,openmars_nlev))
    allocate(openmars_p  (openmars_nlon,openmars_nlat,openmars_nlev))
    allocate(openmars_ps (openmars_nlon,openmars_nlat              ))

    call fiona_start_input('openmars')
    call fiona_input('openmars', 'lon' , openmars_lon)
    call fiona_input('openmars', 'lat' , openmars_lat)
    call fiona_input('openmars', 'lev' , openmars_lev)
    call fiona_input('openmars', 'u'   , openmars_u  )
    call fiona_input('openmars', 'v'   , openmars_v  )
    call fiona_input('openmars', 'temp', openmars_t  )
    call fiona_input('openmars', 'ps'  , openmars_ps )
    call fiona_end_input('openmars')

    allocate(tmp(openmars_nlon,openmars_nlat,openmars_nlev))

    ! Reverse latitude order (from South Pole to North Pole).
    tmp(1,:,1) = openmars_lat(openmars_nlat:1:-1); openmars_lat = tmp(1,:,1)
    tmp(:,:,:) = openmars_u (:,openmars_nlat:1:-1,:); openmars_u  = tmp(:,:,:)
    tmp(:,:,:) = openmars_v (:,openmars_nlat:1:-1,:); openmars_v  = tmp(:,:,:)
    tmp(:,:,:) = openmars_t (:,openmars_nlat:1:-1,:); openmars_t  = tmp(:,:,:)
    tmp(:,:,1) = openmars_ps(:,openmars_nlat:1:-1  ); openmars_ps = tmp(:,:,1)

    ! Flip along the Meridian and change longitude to 0-360.
    tmp(1:openmars_nlon/2,1,1) = openmars_lon(openmars_nlon/2+1:openmars_nlon)
    tmp(openmars_nlon/2+1:openmars_nlon,1,1) = 360 + openmars_lon(1:openmars_nlon/2)
    openmars_lon = tmp(:,1,1)
    tmp(1:openmars_nlon/2,:,:) = openmars_u(openmars_nlon/2+1:openmars_nlon,:,:)
    tmp(openmars_nlon/2+1:openmars_nlon,:,:) = openmars_u(1:openmars_nlon/2,:,:)
    openmars_u = tmp
    tmp(1:openmars_nlon/2,:,:) = openmars_v(openmars_nlon/2+1:openmars_nlon,:,:)
    tmp(openmars_nlon/2+1:openmars_nlon,:,:) = openmars_v(1:openmars_nlon/2,:,:)
    openmars_v = tmp
    tmp(1:openmars_nlon/2,:,:) = openmars_t(openmars_nlon/2+1:openmars_nlon,:,:)
    tmp(openmars_nlon/2+1:openmars_nlon,:,:) = openmars_t(1:openmars_nlon/2,:,:)
    openmars_t = tmp
    tmp(1:openmars_nlon/2,:,1) = openmars_ps(openmars_nlon/2+1:openmars_nlon,:)
    tmp(openmars_nlon/2+1:openmars_nlon,:,1) = openmars_ps(1:openmars_nlon/2,:)
    openmars_ps = tmp(:,:,1)

    ! Reverse vertical levels (from top to bottom).
    tmp = openmars_u(:,:,openmars_nlev:1:-1); openmars_u = tmp
    tmp = openmars_v(:,:,openmars_nlev:1:-1); openmars_v = tmp
    tmp = openmars_t(:,:,openmars_nlev:1:-1); openmars_t = tmp
    tmp(1,1,:) = openmars_lev(openmars_nlev:1:-1); openmars_lev = tmp(1,1,:)

    deallocate(tmp)

    ! Calcuate pressure based on sigma formula.
    do k = 1, openmars_nlev
      openmars_p(:,:,k) = openmars_lev(k) * openmars_ps(:,:)
    end do

  end subroutine openmars_reader_run

  subroutine openmars_reader_final()

    if (allocated(openmars_lon)) deallocate(openmars_lon)
    if (allocated(openmars_lat)) deallocate(openmars_lat)
    if (allocated(openmars_lev)) deallocate(openmars_lev)
    if (allocated(openmars_u  )) deallocate(openmars_u  )
    if (allocated(openmars_v  )) deallocate(openmars_v  )
    if (allocated(openmars_t  )) deallocate(openmars_t  )
    if (allocated(openmars_p  )) deallocate(openmars_p  )
    if (allocated(openmars_ps )) deallocate(openmars_ps )

  end subroutine openmars_reader_final

end module openmars_reader_mod
