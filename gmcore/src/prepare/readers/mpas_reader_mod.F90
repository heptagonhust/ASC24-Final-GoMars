module mpas_reader_mod

  use fiona
  use flogger
  use const_mod
  use formula_mod
  use block_mod
  use process_mod

  implicit none

  ! This is not the raw MPAS data.
  integer mpas_nlon
  integer mpas_nlat
  integer mpas_nlev

  real(r8), allocatable, dimension(:    ) :: mpas_lon
  real(r8), allocatable, dimension(:    ) :: mpas_lat
  real(r8), allocatable, dimension(:,:,:) :: mpas_u
  real(r8), allocatable, dimension(:,:,:) :: mpas_v
  real(r8), allocatable, dimension(:,:,:) :: mpas_pt
  real(r8), allocatable, dimension(:,:,:) :: mpas_t
  real(r8), allocatable, dimension(:,:,:) :: mpas_p
  real(r8), allocatable, dimension(:,:  ) :: mpas_ps
  real(r8), allocatable, dimension(:,:  ) :: mpas_zs

contains

  subroutine mpas_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    real(r8), allocatable :: tmp(:,:,:)
    integer k

    call mpas_reader_final()

    if (proc%is_root()) call log_notice('Use MPAS ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('mpas', file_path=bkg_file)
    call fiona_get_dim('mpas', 'lon', size=mpas_nlon)
    call fiona_get_dim('mpas', 'lat', size=mpas_nlat)
    call fiona_get_dim('mpas', 'lev', size=mpas_nlev)

    allocate(tmp      (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_lon (mpas_nlon))
    allocate(mpas_lat (mpas_nlat))
    allocate(mpas_u   (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_v   (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_pt  (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_t   (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_p   (mpas_nlon,mpas_nlat,mpas_nlev))
    allocate(mpas_ps  (mpas_nlon,mpas_nlat))
    allocate(mpas_zs  (mpas_nlon,mpas_nlat))

    call fiona_start_input('mpas')
    call fiona_input('mpas', 'lon', mpas_lon )
    call fiona_input('mpas', 'lat', mpas_lat )
    call fiona_input('mpas', 'u'  , mpas_u   ); tmp = mpas_u (:,:,mpas_nlev:1:-1); mpas_u  = tmp
    call fiona_input('mpas', 'v'  , mpas_v   ); tmp = mpas_v (:,:,mpas_nlev:1:-1); mpas_v  = tmp
    call fiona_input('mpas', 'pt' , mpas_pt  ); tmp = mpas_pt(:,:,mpas_nlev:1:-1); mpas_pt = tmp
    call fiona_input('mpas', 'ph' , mpas_p   ); tmp = mpas_p (:,:,mpas_nlev:1:-1); mpas_p  = tmp
    call fiona_input('mpas', 'phs', mpas_ps  )
    call fiona_input('mpas', 'ter', mpas_zs  )
    call fiona_end_input('mpas')

    mpas_t = temperature(mpas_pt, mpas_p, 0.0_r8)

    deallocate(tmp)

  end subroutine mpas_reader_run

  subroutine mpas_reader_final()

    if (allocated(mpas_lon )) deallocate(mpas_lon )
    if (allocated(mpas_lat )) deallocate(mpas_lat )
    if (allocated(mpas_u   )) deallocate(mpas_u   )
    if (allocated(mpas_v   )) deallocate(mpas_v   )
    if (allocated(mpas_pt  )) deallocate(mpas_pt  )
    if (allocated(mpas_p   )) deallocate(mpas_p   )
    if (allocated(mpas_ps  )) deallocate(mpas_ps  )
    if (allocated(mpas_zs  )) deallocate(mpas_zs  )

  end subroutine mpas_reader_final

end module mpas_reader_mod
