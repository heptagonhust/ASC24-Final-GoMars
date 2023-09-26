module pole_damp_mod

  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod

  implicit none

  private

  public pole_damp_init
  public pole_damp_final
  public pole_damp_run

  real(r8), allocatable :: c_lon(:)
  real(r8), allocatable :: c_lat(:)

contains

  subroutine pole_damp_init()

    integer j

    call pole_damp_final()

    allocate(c_lon(global_mesh%full_jds:global_mesh%full_jde))
    allocate(c_lat(global_mesh%half_jds:global_mesh%half_jde))

    do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
      c_lon(j) = exp_two_values(pole_damp_coef, 0.0_r8       , &
        real(abs(global_mesh%full_lat_deg(2)), r8)           , &
        real(abs(global_mesh%full_lat_deg(pole_damp_j0)), r8), &
        real(abs(global_mesh%full_lat_deg(j)), r8))
      c_lon(j) = merge(0.0_r8, c_lon(j), c_lon(j) < 1.0e-15)
    end do

    do j = global_mesh%half_jds, global_mesh%half_jde
      c_lat(j) = exp_two_values(pole_damp_coef, 0.0_r8       , &
        real(abs(global_mesh%half_lat_deg(1)), r8)           , &
        real(abs(global_mesh%half_lat_deg(pole_damp_j0)), r8), &
        real(abs(global_mesh%half_lat_deg(j)), r8))
      c_lat(j) = merge(0.0_r8, c_lat(j), c_lat(j) < 1.0e-15)
    end do

  end subroutine pole_damp_init

  subroutine pole_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

  end subroutine pole_damp_final

  subroutine pole_damp_run(block, dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(8), intent(in) :: dt

    integer i, j, k
    real(r8) tmp(global_mesh%full_kms:global_mesh%full_kme)

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, &
               v    => dstate%v_lat)
    ! Set upper and lower boundary conditions.
    do k = mesh%full_kds - 1, mesh%full_kms, -1
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          ! u(i,j,k) = 2 * u(i,j,k+1) - u(i,j,k+2)
          u(i,j,k) = 3 * u(i,j,k+1) - 3 * u(i,j,k+2) + u(i,j,k+3)
        end do
      end do
    end do
    do k = mesh%full_kde + 1, mesh%full_kme
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          ! u(i,j,k) = 2 * u(i,j,k-1) - u(i,j,k-2)
          u(i,j,k) = 3 * u(i,j,k-1) - 3 * u(i,j,k-2) + u(i,j,k-3)
        end do
      end do
    end do
    do k = mesh%full_kds - 1, mesh%full_kms, -1
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          ! v(i,j,k) = 2 * v(i,j,k+1) - v(i,j,k+2)
          v(i,j,k) = 3 * v(i,j,k+1) - 3 * v(i,j,k+2) + v(i,j,k+3)
        end do
      end do
    end do
    do k = mesh%full_kde + 1, mesh%full_kme
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          ! v(i,j,k) = 2 * v(i,j,k-1) - v(i,j,k-2)
          v(i,j,k) = 3 * v(i,j,k-1) - 3 * v(i,j,k-2) + v(i,j,k-3)
        end do
      end do
    end do

    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids, mesh%half_ide
        tmp = u(i,j,:)
        do k = mesh%full_kds, mesh%full_kde
          u(i,j,k) = c_lon(j) * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c_lon(j)) * tmp(k)
        end do
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_jds, mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        tmp = v(i,j,:)
        do k = mesh%full_kds, mesh%full_kde
          v(i,j,k) = c_lat(j) * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c_lat(j)) * tmp(k)
        end do
      end do
    end do
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine pole_damp_run

end module pole_damp_mod
