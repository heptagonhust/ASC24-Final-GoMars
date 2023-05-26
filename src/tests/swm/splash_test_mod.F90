module splash_test_mod

  use const_mod
  use latlon_parallel_mod
  use block_mod

  implicit none

  private

  public splash_test_set_params
  public splash_test_set_ic

  real(r8) gh0
  real(r8) gh1
  real(r8), parameter :: R = 500.0d3

contains

  subroutine splash_test_set_params()

    omega = 0
    gh0 = 50 * g
    gh1 = g

  end subroutine splash_test_set_params

  subroutine splash_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) d
    integer i, j

    associate (mesh   => block%mesh           , &
               u      => block%dstate(1)%u_lon, &
               v      => block%dstate(1)%v_lat, &
               gz     => block%dstate(1)%gz   , &
               gzs    => block%static%gzs)
    u = 0
    v = 0
    gzs = 0

    do j = mesh%full_jds, mesh%full_jde
      d = radius * (pi05 - mesh%full_lat(j))
      do i = mesh%full_ids, mesh%full_ide
        gz(i,j,1) = merge(real(gh0 + gh1 * cos(pi05 * d / R), r8), gh0, d < R)
      end do
    end do
    call fill_halo(block%halo, gz, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine splash_test_set_ic

end module splash_test_mod
