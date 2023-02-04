module mountain_wave_test_mod

  use flogger
  use namelist_mod
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use operators_mod
  use diag_state_mod

  implicit none

  private

  public mountain_wave_test_set_diag
  public mountain_wave_test_set_ic

  real(r8), parameter :: T0   = 288.d0      ! K
  real(r8), parameter :: h0   = 2000.d0     ! m
  real(r8), parameter :: d    = 1.5e6 
  real(r8), parameter :: u0   = 20.d0       ! m s-1
  real(r8), parameter :: lonc = pi05
  real(r8), parameter :: latc = pi / 6.0
  real(r8), parameter :: kap  = 2.d0 / 7.d0
  real(r8), parameter :: psp  = 93000.d0    ! Pa
  real(r8), parameter :: N    = 0.0182      ! s-1

contains

  subroutine mountain_wave_test_set_diag(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call diag_state(iblk)%init_pressure_levels(blocks(iblk), [700.0e2_r8], instance)
    end do

  end subroutine mountain_wave_test_set_diag

  subroutine mountain_wave_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               phs    => block%dstate(1)%phs   , &
               ph_lev => block%dstate(1)%ph_lev, &
               ph     => block%dstate(1)%ph    , &
               t      => block%dstate(1)%t     , &
               pt     => block%dstate(1)%pt    , &
               gzs    => block%static%gzs)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          u(i,j,k) = u0 * cos_lat
        end do
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    v = 0.0

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        full_lon = mesh%full_lon(i)
        r = radius * acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        gzs(i,j) = g * h0 * exp(-(r / d)**2)
      end do
    end do
    call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        phs(i,j) = psp * exp(-0.5_r8 * radius * N**2 * u0 / g**2 / kap * (u0 / radius + 2.0_r8 * omega) * &
                   (sin_lat**2 - 1.0_r8) - N**2 / g**2 / kap * gzs(i,j))
      end do
    end do
    call fill_halo(block%halo, phs, full_lon=.true., full_lat=.true.)

    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
        end do
      end do
    end do
    call fill_halo(block%halo, ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph(i,j,k) = 0.5d0 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
        end do
      end do
    end do
    call fill_halo(block%halo, ph, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          t (i,j,k) = 288.d0
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(block%halo, t, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    if (nonhydrostatic) then
      call calc_gz_lev(block, block%dstate(1))
    end if
    end associate
  
  end subroutine mountain_wave_test_set_ic
  
end module mountain_wave_test_mod
