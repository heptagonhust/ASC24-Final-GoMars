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
  real(r8), allocatable, dimension(:,:,:) :: era5_qv
  real(r8), allocatable, dimension(:,:,:) :: era5_ql
  real(r8), allocatable, dimension(:,:,:) :: era5_qi
  real(r8), allocatable, dimension(:,:,:) :: era5_qr
  real(r8), allocatable, dimension(:,:,:) :: era5_qs
  real(r8), allocatable, dimension(:,:  ) :: era5_ps
  real(r8), allocatable, dimension(:,:  ) :: era5_zs

contains

  subroutine era5_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer i, j, k, k0
    real(r8) qm
    real(r8), allocatable :: tmp(:,:)

    call era5_reader_final()

    if (proc%is_root()) call log_notice('Use ERA5 ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('era5', file_path=bkg_file)
    call fiona_get_dim('era5', 'longitude', size=era5_nlon)
    call fiona_get_dim('era5', 'latitude' , size=era5_nlat)
    call fiona_get_dim('era5', 'level'    , size=era5_nlev)

    allocate(era5_lon(era5_nlon))
    allocate(era5_lat(era5_nlat))
    allocate(era5_lev(era5_nlev))
    allocate(era5_u  (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_v  (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_t  (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_qv (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_ps (era5_nlon,era5_nlat          ))
    allocate(era5_zs (era5_nlon,era5_nlat          ))

    call fiona_start_input('era5')
    call fiona_input('era5', 'longitude', era5_lon)
    call fiona_input('era5', 'latitude' , era5_lat)
    call fiona_input('era5', 'level'    , era5_lev)
    call fiona_input('era5', 'u'        , era5_u  )
    call fiona_input('era5', 'v'        , era5_v  )
    call fiona_input('era5', 't'        , era5_t  )
    call fiona_input('era5', 'q'        , era5_qv )
    call fiona_input('era5', 'sp'       , era5_ps )
    call fiona_input('era5', 'zs'       , era5_zs )
    if (fiona_has_var('era5', 'clwc')) then
      allocate(era5_ql(era5_nlon,era5_nlat,era5_nlev))
      call fiona_input('era5', 'clwc', era5_ql)
    end if
    if (fiona_has_var('era5', 'ciwc')) then
      allocate(era5_qi(era5_nlon,era5_nlat,era5_nlev))
      call fiona_input('era5', 'ciwc', era5_qi)
    end if
    if (fiona_has_var('era5', 'crwc')) then
      allocate(era5_qr(era5_nlon,era5_nlat,era5_nlev))
      call fiona_input('era5', 'crwc', era5_qr)
    end if
    if (fiona_has_var('era5', 'cswc')) then
      allocate(era5_qs(era5_nlon,era5_nlat,era5_nlev))
      call fiona_input('era5', 'cswc', era5_qs)
    end if
    call fiona_end_input('era5')

    ! Reverse latitude order (from South Pole to North Pole).
    allocate(tmp(era5_nlon,era5_nlat))
    tmp(1,:) = era5_lat(era5_nlat:1:-1); era5_lat = tmp(1,:)
    do k = 1, era5_nlev
      tmp = era5_u (:,era5_nlat:1:-1,k); era5_u (:,:,k) = tmp
      tmp = era5_v (:,era5_nlat:1:-1,k); era5_v (:,:,k) = tmp
      tmp = era5_t (:,era5_nlat:1:-1,k); era5_t (:,:,k) = tmp
      tmp = era5_qv(:,era5_nlat:1:-1,k); era5_qv(:,:,k) = tmp
      if (allocated(era5_ql)) then
        tmp = era5_ql(:,era5_nlat:1:-1,k); era5_ql(:,:,k) = tmp
      end if
      tmp = era5_qi(:,era5_nlat:1:-1,k); era5_qi(:,:,k) = tmp
      tmp = era5_qr(:,era5_nlat:1:-1,k); era5_qr(:,:,k) = tmp
      tmp = era5_qs(:,era5_nlat:1:-1,k); era5_qs(:,:,k) = tmp
    end do
    tmp = era5_ps (:,era5_nlat:1:-1  ); era5_ps = tmp
    tmp = era5_zs (:,era5_nlat:1:-1  ); era5_zs = tmp
    deallocate(tmp)

    ! Change units.
    era5_lev = era5_lev * 100.0_r8
    era5_zs  = era5_zs / g

    ! Change moist mixing ratio to dry mixing ratio.
    do k = 1, era5_nlev
      do j = 1, era5_nlat
        do i = 1, era5_nlon
          qm = era5_qv(i,j,k) + era5_ql(i,j,k) + era5_qi(i,j,k) + era5_qr(i,j,k) + era5_qs(i,j,k)
          era5_qv(i,j,k) = era5_qv(i,j,k) / (1 - qm)
          era5_ql(i,j,k) = era5_ql(i,j,k) / (1 - qm)
          era5_qi(i,j,k) = era5_qi(i,j,k) / (1 - qm)
          era5_qr(i,j,k) = era5_qr(i,j,k) / (1 - qm)
          era5_qs(i,j,k) = era5_qs(i,j,k) / (1 - qm)
        end do
      end do
    end do

  end subroutine era5_reader_run

  subroutine era5_reader_final()

    if (allocated(era5_lon)) deallocate(era5_lon)
    if (allocated(era5_lat)) deallocate(era5_lat)
    if (allocated(era5_lev)) deallocate(era5_lev)
    if (allocated(era5_u  )) deallocate(era5_u  )
    if (allocated(era5_v  )) deallocate(era5_v  )
    if (allocated(era5_t  )) deallocate(era5_t  )
    if (allocated(era5_qv )) deallocate(era5_qv )
    if (allocated(era5_ql )) deallocate(era5_ql )
    if (allocated(era5_qi )) deallocate(era5_qi )
    if (allocated(era5_qr )) deallocate(era5_qr )
    if (allocated(era5_qs )) deallocate(era5_qs )
    if (allocated(era5_ps )) deallocate(era5_ps )
    if (allocated(era5_zs )) deallocate(era5_zs )

  end subroutine era5_reader_final

end module era5_reader_mod
