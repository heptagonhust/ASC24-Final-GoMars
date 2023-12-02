module dcmip12_test_mod

  use const_mod, only: r8, pi, pi2, rd, g, radius
  use namelist_mod, dt => dt_adv
  use latlon_parallel_mod
  use history_mod
  use block_mod
  use vert_coord_mod
  use operators_mod
  use adv_mod
  use tracer_mod

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
    rho0 = p0 / (rd * T0)

  end subroutine dcmip12_test_init

  subroutine dcmip12_test_set_ic()

    real(r8) z
    integer i, j, k, iblk

    call tracer_add('dcmip12', dt, 'q0', 'background tracer')
    call tracer_add('dcmip12', dt, 'q1', 'test tracer'      )

    call tracer_allocate()

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 gz    => blocks(iblk)%dstate(1)%gz, &
                 q     => tracers(iblk)%q          )
      ! Background
      q%d(:,:,:,1) = 1
      ! Test
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            z = gz%d(i,j,k) / g
            if (z1 < z .and. z < z2) then
              q%d(i,j,k,2) = 0.5_r8 * (1 + cos(pi2 * (z - z0) / (z2 - z1)))
            else
              q%d(i,j,k,2) = 0
            end if
          end do
        end do
      end do
      call fill_halo(q, 2, cross_pole=.true.)
      end associate
    end do

  end subroutine dcmip12_test_set_ic

  subroutine dcmip12_test_set_uv(time_in_seconds, itime)

    real(r8), intent(in) :: time_in_seconds
    integer, intent(in) :: itime

    integer iblk, i, j, k
    real(r8) lon, lat, rho, cos_t, dmgdlev

    cos_t = cos(pi * time_in_seconds / tau)

    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                      , &
                 dstate  => blocks(iblk)%dstate(itime)        , &
                 mesh    => blocks(iblk)%mesh                 , &
                 mgs     => blocks(iblk)%dstate(itime)%mgs    , &
                 mg_lev  => blocks(iblk)%dstate(itime)%mg_lev , &
                 mg      => blocks(iblk)%dstate(itime)%mg     , &
                 dmg_lon => blocks(iblk)%aux%dmg_lon          , &
                 dmg_lat => blocks(iblk)%aux%dmg_lat          , &
                 dmg_lev => blocks(iblk)%dstate(itime)%dmg_lev, &
                 ph_lev  => blocks(iblk)%dstate(itime)%ph_lev , &
                 ph      => blocks(iblk)%dstate(itime)%ph     , &
                 tv      => blocks(iblk)%dstate(itime)%tv     , &
                 gz_lev  => blocks(iblk)%dstate(itime)%gz_lev , &
                 gz      => blocks(iblk)%dstate(itime)%gz     , &
                 u       => blocks(iblk)%dstate(itime)%u_lon  , &
                 v       => blocks(iblk)%dstate(itime)%v_lat  , &
                 mfx_lon => blocks(iblk)%aux%mfx_lon          , &
                 mfy_lat => blocks(iblk)%aux%mfy_lat          , &
                 we      => blocks(iblk)%dstate(itime)%we_lev )
      mgs%d = p0
      call calc_mg(block, dstate)
      ph_lev%d = mg_lev%d
      ph%d = mg%d
      tv%d = T0
      call calc_dmg(block, dstate)
      call calc_gz_lev(block, dstate)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          lat = mesh%full_lat(j)
          do i = mesh%half_ids, mesh%half_ide
            lon = mesh%half_lon(i)
            u%d(i,j,k) = u0 * cos(lat)
          end do
        end do
      end do
      call fill_halo(u)
      mfx_lon%d = u%d * dmg_lon%d
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          lat = mesh%half_lat(j)
          do i = mesh%full_ids, mesh%full_ide
            lon = mesh%full_lon(i)
            rho = mg%d(i,j,k) / (rd * T0)
            v%d(i,j,k) = -radius * w0 * pi * rho0 / (K0 * ztop * rho) * &
              cos(lat) * sin(K0 * lat) * cos(pi * gz%d(i,j,k) / (g * ztop)) * cos_t
          end do
        end do
      end do
      call fill_halo(v)
      mfy_lat%d = v%d * dmg_lat%d
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          lat = mesh%full_lat(j)
          do i = mesh%full_ids, mesh%full_ide
            lon = mesh%full_lon(i)
            rho = mg_lev%d(i,j,k) / (rd * T0)
            dmgdlev = dmg_lev%d(i,j,k) / mesh%half_dlev(k)
            we%d(i,j,k) = -dmgdlev * g * w0 * rho0 / K0 * (                   &
              -2 * sin(K0 * lat) * sin(lat) + K0 * cos(lat) * cos(K0 * lat) &
            ) * sin(pi * gz_lev%d(i,j,k) / (g * ztop)) * cos_t / p0
          end do
        end do
      end do
      end associate
    end do

  end subroutine dcmip12_test_set_uv

end module dcmip12_test_mod
