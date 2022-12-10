module dcmip12_test_mod

  use const_mod, only: r8, pi, pi2, Rd, g, radius
  use namelist_mod, dt => dt_adv
  use parallel_mod
  use history_mod
  use block_mod
  use vert_coord_mod
  use operators_mod
  use interp_mod
  use adv_mod

  implicit none

  private

  public dcmip12_test_init
  public dcmip12_test_set_ic
  public dcmip12_test_set_uv

  real(r8), parameter :: T0   = 300
  real(r8), parameter :: p0   = 1.0e5_r8
  real(r8), parameter :: K0   = 5
  real(r8), parameter :: u0   = 40
  real(r8), parameter :: w0   = 0.15
  real(r8), parameter :: z1   = 2000
  real(r8), parameter :: z2   = 5000
  real(r8), parameter :: z0   = 0.5_r8 * (z1 + z2)
  real(r8), parameter :: ztop = 12000
  real(r8), parameter :: tau  = 86400

  real(r8) rho0

contains

  subroutine dcmip12_test_init()

    ptop = 25494.4
    rho0 = p0 / (Rd * T0)

  end subroutine dcmip12_test_init

  subroutine dcmip12_test_set_ic(block)

    type(block_type), intent(inout) :: block

    real(r8) z
    integer i, j, k

    call adv_add_tracer('dcmip12', dt, 'q0', 'background tracer')
    call adv_add_tracer('dcmip12', dt, 'q1', 'test tracer'      )

    call adv_allocate_tracers(block)

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               gz   => block%dstate(1)%gz      , &
               q    => block%adv_batches(1)%q  )
    ! Background
    q(:,:,:,1,old) = 1
    ! Test
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          z = gz(i,j,k) / g
          if (z1 < z .and. z < z2) then
            q(i,j,k,2,old) = 0.5_r8 * (1 + cos(pi2 * (z - z0) / (z2 - z1)))
          else
            q(i,j,k,2,old) = 0
          end if
        end do
      end do
    end do
    call fill_halo(block%halo, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine dcmip12_test_set_ic

  subroutine dcmip12_test_set_uv(block, dstate, time_in_seconds)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: time_in_seconds

    integer i, j, k
    real(r8) lon, lat, rho, cos_t, dphdlev

    cos_t = cos(pi * time_in_seconds / tau)

    associate (mesh    => block%mesh    , &
               phs     => dstate%phs    , &
               ph_lev  => dstate%ph_lev , &
               ph      => dstate%ph     , &
               m       => dstate%m      , &
               m_lon   => dstate%m_lon  , &
               m_lat   => dstate%m_lat  , &
               m_lev   => dstate%m_lev  , &
               t       => dstate%t      , &
               gz_lev  => dstate%gz_lev , &
               gz      => dstate%gz     , &
               u       => dstate%u_lon  , &
               v       => dstate%v_lat  , &
               mfx_lon => dstate%mfx_lon, &
               mfy_lat => dstate%mfy_lat, &
               we      => dstate%we_lev)
    phs = p0
    t   = T0
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph(i,j,k) = 0.5d0 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
        end do
      end do
    end do
    call calc_m(block, dstate)
    call calc_gz_lev(block, dstate)
    call interp_lev_edge_to_cell(mesh, gz_lev, gz)
    ! Set invariant pressure thickness.
    do k = mesh%full_kds, mesh%full_kde
      m    (:,:,k) = p0 * mesh%full_dlev(k)
      m_lon(:,:,k) = p0 * mesh%full_dlev(k)
      m_lat(:,:,k) = p0 * mesh%full_dlev(k)
    end do
    do k = mesh%half_kds, mesh%half_kde
      m_lev(:,:,k) = p0 * mesh%half_dlev(k)
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        lat = mesh%full_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          lon = mesh%half_lon(i)
          u(i,j,k) = u0 * cos(lat)
          mfx_lon(i,j,k) = u(i,j,k) * m_lon(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        lat = mesh%half_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          rho = ph(i,j,k) / (Rd * T0)
          v(i,j,k) = -radius * w0 * pi * rho0 / (K0 * ztop * rho) * &
            cos(lat) * sin(K0 * lat) * cos(pi * gz(i,j,k) / (g * ztop)) * cos_t
          mfy_lat(i,j,k) = v(i,j,k) * m_lat(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          rho = ph_lev(i,j,k) / (Rd * T0)
          dphdlev = m_lev(i,j,k) / mesh%half_dlev(k)
          we(i,j,k) = -dphdlev * g * w0 * rho0 / K0 * (                   &
            -2 * sin(K0 * lat) * sin(lat) + K0 * cos(lat) * cos(K0 * lat) &
          ) * sin(pi * gz_lev(i,j,k) / (g * ztop)) * cos_t / p0
        end do
      end do
    end do
    end associate

  end subroutine dcmip12_test_set_uv

end module dcmip12_test_mod
