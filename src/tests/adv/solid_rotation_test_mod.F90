module solid_rotation_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use latlon_parallel_mod
  use block_mod
  use adv_mod
  use tracer_mod

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

  subroutine solid_rotation_test_set_ic()

    integer iblk, i, j
    real(8) lon, lat, r, R0

    call tracer_add('solid_rotation', dt, 'q0', 'background tracer' )
    call tracer_add('solid_rotation', dt, 'q1', 'cosine bell tracer')

    call tracer_allocate()

    R0 = radius / 3.0_r8

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)     , &
                 mesh  => blocks(iblk)%mesh, &
                 q     => tracers(iblk)%q  )
      ! Background
      q%d(:,:,:,1) = 1
      ! Cosine bell
      do j = mesh%full_jds, mesh%full_jde
        lat = mesh%full_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          r = great_circle(radius, lon0, lat0, lon, lat)
          if (r < R0) then
            q%d(i,j,1,2) = h0 / 2.0_r8 * (1 + cos(pi * r / R0))
          else
            q%d(i,j,1,2) = 0
          end if
        end do
      end do
      call fill_halo(q, 2)
      end associate
    end do

  end subroutine solid_rotation_test_set_ic

  subroutine solid_rotation_test_set_uv(time_in_seconds, itime)

    real(r8), intent(in) :: time_in_seconds
    integer, intent(in) :: itime

    integer iblk, i, j
    real(r8) lon, lat

    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                    , &
                 mesh    => blocks(iblk)%mesh               , &
                 dmg     => blocks(iblk)%dstate(itime)%dmg  , &
                 dmg_lon => blocks(iblk)%aux%dmg_lon        , &
                 dmg_lat => blocks(iblk)%aux%dmg_lat        , &
                 u       => blocks(iblk)%dstate(itime)%u_lon, &
                 v       => blocks(iblk)%dstate(itime)%v_lat, &
                 mfx     => blocks(iblk)%aux%mfx_lon        , &
                 mfy     => blocks(iblk)%aux%mfy_lat        )
      dmg%d = 1; dmg_lon%d = 1; dmg_lat%d = 1
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        lat = mesh%full_lat(j)
        do i = mesh%half_ids, mesh%half_ide
          lon = mesh%half_lon(i)
          u%d(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
        end do
      end do
      call fill_halo(u)
      mfx%d = u%d
      do j = mesh%half_jds, mesh%half_jde
        lat = mesh%half_lat(j)
        do i = mesh%full_ids, mesh%full_ide
          lon = mesh%full_lon(i)
          v%d(i,j,1) = -u0 * sin(lon) * sin(alpha)
        end do
      end do
      call fill_halo(v)
      mfy%d = v%d
      end associate
    end do

  end subroutine solid_rotation_test_set_uv

end module solid_rotation_test_mod
