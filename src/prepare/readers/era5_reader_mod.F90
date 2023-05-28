module era5_reader_mod

  use fiona
  use flogger
  use const_mod
  use namelist_mod
  use formula_mod
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
  real(r8), allocatable, dimension(:,:,:) :: era5_z
  real(r8), allocatable, dimension(:,:,:) :: era5_qv
  real(r8), allocatable, dimension(:,:,:) :: era5_ql
  real(r8), allocatable, dimension(:,:,:) :: era5_qi
  real(r8), allocatable, dimension(:,:,:) :: era5_qr
  real(r8), allocatable, dimension(:,:,:) :: era5_qs
  real(r8), allocatable, dimension(:,:,:) :: era5_pd
  real(r8), allocatable, dimension(:,:  ) :: era5_ps
  real(r8), allocatable, dimension(:,:  ) :: era5_psd
  real(r8), allocatable, dimension(:,:  ) :: era5_zs

contains

  subroutine era5_reader_run(min_lon, max_lon, min_lat, max_lat)

    use mpi

    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat

    integer i, j, k, k0
    real(r8) qm1, qm2, qm, tv1, tv2, sum_q_lev, dz, rho

    call era5_reader_final()

    if (proc%is_root()) call log_notice('Use ERA5 ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('era5', file_path=bkg_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('era5', 'longitude', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('era5', 'latitude', span=[90, -90], flip=.true.)
    call fiona_get_dim('era5', 'level', size=era5_nlev)
    allocate(era5_lev(era5_nlev))
    call fiona_start_input('era5')
    call fiona_input_range('era5', 'longitude', era5_lon, coord_range=[min_lon, max_lon]); era5_nlon = size(era5_lon)
    call fiona_input_range('era5', 'latitude' , era5_lat, coord_range=[min_lat, max_lat]); era5_nlat = size(era5_lat)
    call fiona_input      ('era5', 'level'    , era5_lev); era5_nlev = size(era5_lev)
    call fiona_input_range('era5', 'u'        , era5_u  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'v'        , era5_v  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 't'        , era5_t  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'z'        , era5_z  , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'q'        , era5_qv , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'sp'       , era5_ps , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_input_range('era5', 'zs'       , era5_zs , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'clwc')) &
    call fiona_input_range('era5', 'clwc'     , era5_ql , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'ciwc')) &
    call fiona_input_range('era5', 'ciwc'     , era5_qi , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'crwc')) &
    call fiona_input_range('era5', 'crwc'     , era5_qr , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
  if (fiona_has_var('era5', 'cswc')) &
    call fiona_input_range('era5', 'cswc'     , era5_qs , coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_end_input('era5')

    ! Change units.
    era5_lev = era5_lev * 100.0_r8
    era5_z   = era5_z  / g
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

    if (proc%is_root()) call log_notice('Calculate ERA5 dry-air pressure or weight.')
    allocate(era5_pd (era5_nlon,era5_nlat,era5_nlev))
    allocate(era5_psd(era5_nlon,era5_nlat))
    do j = 1, era5_nlat
      do i = 1, era5_nlon
        k0 = 1
        do k = 1, era5_nlev
          if (era5_lev(k) >= era5_ps(i,j)) then
            k0 = k - 1
            exit
          end if
        end do
        k0 = min(k0, era5_nlev - 1)
        era5_pd(i,j,:) = era5_lev(:)
        sum_q_lev = 0
        ! Pressure levels
        do k = 2, k0
          qm1 = era5_qv(i,j,k-1) + era5_ql(i,j,k-1) + era5_qi(i,j,k-1) + era5_qr(i,j,k-1) + era5_qs(i,j,k-1)
          qm2 = era5_qv(i,j,k  ) + era5_ql(i,j,k  ) + era5_qi(i,j,k  ) + era5_qr(i,j,k  ) + era5_qs(i,j,k  )
          qm  = 0.5_r8 * (qm1 + qm2)
          tv1 = virtual_temperature(era5_t(i,j,k-1), era5_qv(i,j,k-1), qm1)
          tv2 = virtual_temperature(era5_t(i,j,k  ), era5_qv(i,j,k  ), qm2)
          rho = 0.5_r8 * (era5_lev(k-1) / rd / tv1 + era5_lev(k) / rd / tv2)
          dz = era5_z(i,j,k-1) - era5_z(i,j,k)
          sum_q_lev = sum_q_lev + g * rho * qm / (1 + qm) * dz
          era5_pd(i,j,k) = era5_pd(i,j,k) - sum_q_lev
        end do
        ! Surface
        qm1 = era5_qv(i,j,k0) + era5_ql(i,j,k0) + era5_qi(i,j,k0) + era5_qr(i,j,k0) + era5_qs(i,j,k0)
        qm2 = era5_qv(i,j,k0) + era5_ql(i,j,k0) + era5_qi(i,j,k0) + era5_qr(i,j,k0) + era5_qs(i,j,k0)
        qm  = 0.5_r8 * (qm1 + qm2)
        tv1 = virtual_temperature(era5_t(i,j,k0), era5_qv(i,j,k0), qm1)
        tv2 = virtual_temperature(era5_t(i,j,k0), era5_qv(i,j,k0), qm2)
        rho = 0.5_r8 * (era5_lev(k0) / rd / tv1 + era5_ps(i,j) / rd / tv2)
        dz = era5_z(i,j,k0) - era5_zs(i,j)
        sum_q_lev = sum_q_lev + g * rho * qm / (1 + qm) * dz
        era5_psd(i,j) = era5_ps(i,j) - sum_q_lev
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
    if (allocated(era5_z  )) deallocate(era5_z  )
    if (allocated(era5_qv )) deallocate(era5_qv )
    if (allocated(era5_ql )) deallocate(era5_ql )
    if (allocated(era5_qi )) deallocate(era5_qi )
    if (allocated(era5_qr )) deallocate(era5_qr )
    if (allocated(era5_qs )) deallocate(era5_qs )
    if (allocated(era5_pd )) deallocate(era5_pd )
    if (allocated(era5_ps )) deallocate(era5_ps )
    if (allocated(era5_psd)) deallocate(era5_psd)
    if (allocated(era5_zs )) deallocate(era5_zs )

  end subroutine era5_reader_final

end module era5_reader_mod
