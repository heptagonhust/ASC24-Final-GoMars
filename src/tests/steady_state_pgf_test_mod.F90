module steady_state_pgf_test_mod

  use flogger
  use const_mod, only: r8, pi, Rd, g, omega
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public steady_state_pgf_test_set_params
  public steady_state_pgf_test_set_ic

  real(r8), parameter :: T0    = 300.d0      ! K
  real(r8), parameter :: h0    = 2000.d0     ! m
  real(r8), parameter :: p0    = 1.0e5       ! pa
  real(r8), parameter :: lonc  = 3.d0 * pi / 2
  real(r8), parameter :: latc  = 0.0
  real(r8), parameter :: Rm    = 3.d0 * pi / 4
  real(r8), parameter :: gamma = 0.0065d0
  real(r8), parameter :: osm   = pi / 16.d0

contains

  subroutine steady_state_pgf_test_set_params()

    omega = 0.0

  end subroutine steady_state_pgf_test_set_params

  subroutine steady_state_pgf_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r, height
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
    u = 0.0
    v = 0.0

    do j = mesh%full_jds, mesh%full_jde
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        full_lon = mesh%full_lon(i)
        r = acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        if (r < Rm) gzs(i,j) = g * h0 / 2.d0 * (1.d0 + cos(pi * r / Rm)) * cos(pi * r / osm)**2
      end do
    end do
    call fill_halo(block%filter_halo, gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        phs(i,j) = p0 * (1.d0 - gamma / T0 * gzs(i,j) / g)**(g / Rd / gamma) 
      end do
    end do
    call fill_halo(block%halo, phs, full_lon=.true., full_lat=.true.)

    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph_lev(i,j,k) = vert_coord_calc_mg_lev(k, phs(i,j))
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
          t (i,j,k) = T0 * (ph(i,j,k) / p0)**(Rd * gamma / g)
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(block%halo, t, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate
  
  end subroutine steady_state_pgf_test_set_ic

end module steady_state_pgf_test_mod
