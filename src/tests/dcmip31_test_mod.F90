module dcmip31_test_mod

  use flogger
  use namelist_mod
  use const_mod, only: r8, pi, Rd, Rd_o_cpd, cpd, g, radius, omega
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use operators_mod
  use diag_state_mod

  implicit none

  private

  public dcmip31_test_set_params
  public dcmip31_test_set_ic

  real(r8), parameter :: X    = 125.0_r8
  real(r8), parameter :: u0   = 20.0_r8      ! m s-1
  real(r8), parameter :: teq  = 300.0_r8     ! K
  real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
  real(r8), parameter :: d    = 5000.0_r8    ! m
  real(r8), parameter :: lonc = 2 * pi / 3
  real(r8), parameter :: latc = 0.0_r8
  real(r8), parameter :: dpt  = 1.0_r8       ! K
  real(r8), parameter :: lz   = 20000.0_r8   ! m
  real(r8), parameter :: N    = 0.01_r8      ! s-1
  real(r8), parameter :: N2   = N**2
  real(r8)            :: t0
  real(r8), parameter :: p0   = 1.0e5_r8     ! Pa

contains

  subroutine dcmip31_test_set_params()

    omega = 0.0_r8
    radius = radius / X
    t0 = g**2 / N2 / cpd

  end subroutine dcmip31_test_set_params

  subroutine dcmip31_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) ts, cos_2lat, r, local_z, local_ztop

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               w      => block%dstate(1)%w_lev , &
               phs    => block%dstate(1)%phs   , &
               ph_lev => block%dstate(1)%ph_lev, &
               ph     => block%dstate(1)%ph    , &
               pt     => block%dstate(1)%pt    , &
               t      => block%dstate(1)%t     , &
               gz_lev => block%dstate(1)%gz_lev, &
               gz     => block%dstate(1)%gz    , &
               gzs    => block%static%gzs)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            u(i,j,k) = u0 * mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)

      v = 0.0_r8
      gzs = 0.0_r8

      do j = mesh%full_jds, mesh%full_jde
        cos_2lat = cos(2 * mesh%full_lat(j))
        ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
        do i = mesh%full_ids, mesh%full_ide
          phs(i,j) = peq * exp(u0 / (4 * t0 * Rd) * (u0 + 2 * omega * radius) * (cos_2lat - 1)) * &
                     (ts / teq)**(1 / Rd_o_cpd)
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

      if (nonhydrostatic) then
        w = 0.0_r8
        do k = mesh%half_kds, mesh%half_kde
          do j = mesh%full_jds, mesh%full_jde
            cos_2lat = cos(2 * mesh%full_lat(j))
            ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
            do i = mesh%full_ids, mesh%full_ide
              gz_lev(i,j,k) = - g**2 / N2 * log(ts / t0 * ((ph_lev(i,j,k) / phs(i,j))**Rd_o_cpd - 1) + 1)
            end do
          end do
        end do
        call fill_halo(block%halo, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          cos_2lat = cos(2 * mesh%full_lat(j))
          ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
          do i = mesh%full_ids, mesh%full_ide
            pt(i,j,k) = ts * (p0 / phs(i,j))**Rd_o_cpd / (ts / t0 * ((ph(i,j,k) / phs(i,j))**Rd_o_cpd - 1) + 1)
          end do
        end do
      end do

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ! Perturbation
            local_z = 0.5_r8 * (gz_lev(i,j,k+1) + gz_lev(i,j,k)) / g
            local_ztop = gz_lev(i,j,mesh%half_kds) / g 
            r = radius * acos(sin(latc) * mesh%full_sin_lat(j) + cos(latc) * mesh%full_cos_lat(j) * cos(mesh%full_lon(i) - lonc))
            pt(i,j,k) = pt(i,j,k) + dpt * d**2 / (d**2 + r**2) * sin(pi * local_z / local_ztop)
          end do
        end do
      end do
      call fill_halo(block%halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate
  
  end subroutine dcmip31_test_set_ic

end module dcmip31_test_mod
