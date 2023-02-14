module solid_rotation_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use block_mod
  use adv_mod

  implicit none

  public solid_rotation_test_init
  public solid_rotation_test_set_ic
  public solid_rotation_test_set_uv

  real(8), parameter :: period = 12 * 86400
  real(8), parameter :: h0     = 1000 ! m
  real(8), parameter :: lon0   = 3 * pi / 2.0_r8
  real(8), parameter :: lat0   = 0
  real(8), parameter :: alpha  = 90.0_r8 * rad
  real(8) u0

contains

  subroutine solid_rotation_test_init()

    u0 = pi2 * radius / period

  end subroutine solid_rotation_test_init

  subroutine solid_rotation_test_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(8) lon, lat, r, R0

    call adv_add_tracer('solid_rotation', dt, 'q0', 'background tracer')
    call adv_add_tracer('solid_rotation', dt, 'q1', 'cosine bell tracer')

    call adv_allocate_tracers(block)

    R0 = radius / 3.0_r8

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               q    => block%adv_batches(1)%q  )
    ! Background
    q(:,:,:,1,old) = 1
    ! Cosine bell
    do j = mesh%full_jds, mesh%full_jde
      lat = mesh%full_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        lon = mesh%full_lon(i)
        r = calc_distance(lon0, lat0, lon, lat)
        if (r < R0) then
          q(i,j,1,2,old) = h0 / 2.0_r8 * (1 + cos(pi * r / R0))
        else
          q(i,j,1,2,old) = 0
        end if
      end do
    end do
    call fill_halo(block%filter_halo, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine solid_rotation_test_set_ic

  subroutine solid_rotation_test_set_uv(block, dstate, time_in_seconds)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat

    associate (mesh => block%mesh    , &
               dmg  => dstate%dmg    , &
               u    => dstate%u_lon  , &
               v    => dstate%v_lat  , &
               mfx  => dstate%mfx_lon, &
               mfy  => dstate%mfy_lat)
    dmg = 1
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_ids, mesh%half_ide
        lon = mesh%half_lon(i)
        u(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    mfx = u
    do j = mesh%half_jds, mesh%half_jde
      lat = mesh%half_lat(j)
      do i = mesh%full_ids, mesh%full_ide
        lon = mesh%full_lon(i)
        v(i,j,1) = -u0 * sin(lon) * sin(alpha)
      end do
    end do
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    mfy = v
    end associate

  end subroutine solid_rotation_test_set_uv

end module solid_rotation_test_mod
