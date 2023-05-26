module moving_vortices_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use latlon_parallel_mod
  use history_mod
  use block_mod
  use adv_mod
  use tracer_mod

  implicit none

  private

  public moving_vortices_test_init
  public moving_vortices_test_set_ic
  public moving_vortices_test_set_uv

  real(8), parameter :: period = 12 * 86400    ! Rotation period
  real(8), parameter :: rho0   = 3             !
  real(8), parameter :: gamma  = 5             ! Field stiffness parameter
  real(8), parameter :: alpha  = pi05          ! Angle between Rotation axis and the equator
  real(8), parameter :: lonp0  = pi            ! Rotation north pole longitude
  real(8), parameter :: latp0  = pi05 - alpha  ! Rotation north pole latitude
  real(8), parameter :: lonv0  = pi * 1.5_r8   ! Initial vortex longitude
  real(8), parameter :: latv0  = 0             ! Initial vortex latitude
  real(8) u0                                   ! Velocity amplitude
  real(8) lonvr0                               ! Initial vortex longitude in rotated coordinate
  real(8) lonvr                                ! Vortex longitude in rotated coordinate
  real(8) latvr                                ! Vortex latitude  in rotated coordinate
  real(8) lonv                                 ! Vortex longitude
  real(8) latv                                 ! Vortex latitude

contains

  subroutine moving_vortices_test_init()

    u0 = pi2 * radius / (12.0_r8 * 86400.0_r8)

    call rotate(lonp0, latp0, lonv0, latv0, lonvr, latvr)
    lonv = lonv0
    latv = latv0
    lonvr0 = lonvr

  end subroutine moving_vortices_test_init

  subroutine moving_vortices_test_set_ic()

    integer iblk, i, j
    real(8) lon, lat, lonr, latr

    call tracer_add('moving_vortices', dt, 'q0', 'background tracer')
    call tracer_add('moving_vortices', dt, 'q1', 'vortex tracer'    )

    call tracer_allocate()

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk), mesh => blocks(iblk)%mesh)
      ! Background tracer
      tracers(iblk)%q(:,:,:,1) = 1
      ! Vortex tracer
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          call rotate(lonv0, latv0, lon, lat, lonr, latr)
          tracers(iblk)%q(i,j,1,2) = 1 - tanh(rho(latr) / gamma * sin(lonr))
        end do
      end do
      call fill_halo(block%filter_halo, tracers(iblk)%q(:,:,:,2), full_lon=.true., full_lat=.true., full_lev=.true.)
      end associate
    end do

  end subroutine moving_vortices_test_set_ic

  subroutine moving_vortices_test_set_uv(time_in_seconds, itime)

    real(r8), intent(in) :: time_in_seconds
    integer, intent(in) :: itime

    integer iblk, i, j
    real(8) lon, lat, dlon, latr

    lonvr = lonvr0 + u0 / radius * time_in_seconds
    if (lonvr > pi2) lonvr = lonvr - pi2
    call rotate_back(lonp0, latp0, lonv, latv, lonvr, latvr)

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)                      , &
                 mesh  => blocks(iblk)%mesh                 , &
                 dmg   => blocks(iblk)%dstate(itime)%dmg    , &
                 u     => blocks(iblk)%dstate(itime)%u_lon  , &
                 v     => blocks(iblk)%dstate(itime)%v_lat  , &
                 mfx   => blocks(iblk)%dstate(itime)%mfx_lon, &
                 mfy   => blocks(iblk)%dstate(itime)%mfy_lat)
      dmg = 1
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        lat = mesh%full_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          lon = mesh%half_lon(i)
          dlon = lon - lonv
          call rotate(lonv, latv, lon, lat, lat_r=latr)
          u(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha)) + &
                     a_omg(latr) * (sin(latv) * cos(lat) - cos(latv) * cos(dlon) * sin(lat))
        end do
      end do
      call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
      mfx = u
      do j = mesh%half_jds, mesh%half_jde
        lat = mesh%half_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          dlon = lon - lonv
          call rotate(lonv, latv, lon, lat, lat_r=latr)
          v(i,j,1) = -u0 * sin(lon) * sin(alpha) + a_omg(latr) * cos(latv) * sin(dlon)
        end do
      end do
      call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
      mfy = v
      end associate
    end do

  end subroutine moving_vortices_test_set_uv

  real(r8) function rho(lat) result(res)

    real(8), intent(in) :: lat

    res = rho0 * cos(lat)

  end function rho

  real(r8) function a_omg(latr) result(res)

    real(8), intent(in) :: latr

    real(r8), parameter :: c = 1.5_r8 * sqrt(3.0_r8)
    real(r8) r

    r = rho(latr)
    if (abs(r) < eps) then
      res = 0
    else
      res = u0 * c / cosh(r)**2 * tanh(r) / r
    end if

  end function a_omg

end module moving_vortices_test_mod
