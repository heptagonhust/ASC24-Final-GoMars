module era5_reader_mod

  use fiona
  use flogger
  use const_mod
  use process_mod

  implicit none

  integer era5_nlon
  integer era5_nlat
  integer era5_nlev

  real(r8), allocatable, dimension(:    ) :: era5_lon
  real(r8), allocatable, dimension(:    ) :: era5_lat
  real(r8), allocatable, dimension(:    ) :: era5_lev
  real(r8), allocatable, dimension(:,:,:) :: era5_u
  real(r8), allocatable, dimension(:,:,:) :: era5_v
  real(r8), allocatable, dimension(:,:,:) :: era5_t
  real(r8), allocatable, dimension(:,:,:) :: era5_q
  real(r8), allocatable, dimension(:,:  ) :: era5_ps
  real(r8), allocatable, dimension(:,:  ) :: era5_mslp
  real(r8), allocatable, dimension(:,:  ) :: era5_psd  ! Dry surface pressure

contains

  subroutine era5_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer k
    real(r8), allocatable :: tmp(:,:)

    call era5_reader_final()

    if (proc%is_root()) call log_notice('Use ERA5 ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('era5', file_path=bkg_file)
    call fiona_get_dim('era5', 'longitude', size=era5_nlon)
    call fiona_get_dim('era5', 'latitude' , size=era5_nlat)
    call fiona_get_dim('era5', 'level'    , size=era5_nlev)

    allocate(era5_lon (era5_nlon))
    allocate(era5_lat (era5_nlat))
    allocate(era5_lev (era5_nlev))
    allocate(era5_u   (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_v   (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_t   (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_q   (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_ps  (era5_nlon,era5_nlat             ))
    allocate(era5_mslp(era5_nlon,era5_nlat             ))
    allocate(era5_psd (era5_nlon,era5_nlat             ))

    call fiona_start_input('era5')
    call fiona_input    ('era5', 'longitude', era5_lon )
    call fiona_input    ('era5', 'latitude' , era5_lat )
    call fiona_input    ('era5', 'level'    , era5_lev )
    call fiona_input    ('era5', 'u'        , era5_u   )
    call fiona_input    ('era5', 'v'        , era5_v   )
    call fiona_input    ('era5', 't'        , era5_t   )
    if (fiona_has_var('era5', 'q')) then
      call fiona_input  ('era5', 'q'        , era5_q   )
    end if
    call fiona_input    ('era5', 'msl'      , era5_mslp)
    call fiona_end_input('era5')

    ! Reverse latitude order (from South Pole to North Pole).
    allocate(tmp(era5_nlon,era5_nlat))
    tmp(1,:) = era5_lat(era5_nlat:1:-1); era5_lat = tmp(1,:)
    do k = 1, era5_nlev
      tmp = era5_u(:,era5_nlat:1:-1,k); era5_u(:,:,k) = tmp
      tmp = era5_v(:,era5_nlat:1:-1,k); era5_v(:,:,k) = tmp
      tmp = era5_t(:,era5_nlat:1:-1,k); era5_t(:,:,k) = tmp
      tmp = era5_q(:,era5_nlat:1:-1,k); era5_q(:,:,k) = tmp
    end do
    tmp = era5_ps  (:,era5_nlat:1:-1); era5_ps   = tmp
    tmp = era5_mslp(:,era5_nlat:1:-1); era5_mslp = tmp
    deallocate(tmp)

    ! Change units.
    era5_lev = era5_lev * 100.0_r8

    ! Change specific humidity to mixing ratio.
    era5_q = era5_q / (1 - era5_q)

  end subroutine era5_reader_run

  subroutine era5_reader_final()

    if (allocated(era5_lon )) deallocate(era5_lon )
    if (allocated(era5_lat )) deallocate(era5_lat )
    if (allocated(era5_lev )) deallocate(era5_lev )
    if (allocated(era5_u   )) deallocate(era5_u   )
    if (allocated(era5_v   )) deallocate(era5_v   )
    if (allocated(era5_t   )) deallocate(era5_t   )
    if (allocated(era5_q   )) deallocate(era5_q   )
    if (allocated(era5_ps  )) deallocate(era5_ps  )
    if (allocated(era5_mslp)) deallocate(era5_mslp)
    if (allocated(era5_psd )) deallocate(era5_psd )

  end subroutine era5_reader_final

end module era5_reader_mod
