! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module vor_damp_mod

  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use operators_mod
  use filter_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: cx(:,:), cy(:,:)

contains

  subroutine vor_damp_init()

    real(r8) r, lat0
    integer j, k

    call vor_damp_final()

    allocate(cx(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(cy(global_mesh%half_nlat,global_mesh%full_nlev))

    select case (vor_damp_order)
    case (2)
      r = 1
      lat0 = abs(global_mesh%full_lat_deg(2))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
          cx(j,k) = vor_damp_coef2 * global_mesh%full_cos_lat(j)**(r - 1) * &
            exp_two_values(vor_damp_top, 1.0_r8, 1.0_r8, real(vor_damp_k0, r8), real(k, r8)) * &
            exp_two_values(vor_damp_pole_x, 1.0_r8, lat0, vor_damp_lat0, abs(global_mesh%full_lat_deg(j))) * &
            global_mesh%le_lon(j) * global_mesh%de_lon(j) / dt_dyn
        end do
      end do
      lat0 = abs(global_mesh%half_lat_deg(1))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%half_jds, global_mesh%half_jde
          cy(j,k) = vor_damp_coef2 * global_mesh%half_cos_lat(j)**(r - 1) * &
            exp_two_values(vor_damp_top, 1.0_r8, 1.0_r8, real(vor_damp_k0, r8), real(k, r8)) * &
            exp_two_values(vor_damp_pole_y, 1.0_r8, lat0, vor_damp_lat0, abs(global_mesh%half_lat_deg(j))) * &
            global_mesh%le_lat(j) * global_mesh%de_lat(j) / dt_dyn
        end do
      end do
    end select

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(cx)) deallocate(cx)
    if (allocated(cy)) deallocate(cy)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(8), intent(in) :: dt

    integer i, j, k

    call calc_vor(block, dstate)

    associate (mesh => block%mesh       , &
               vor  => block%aux%vor    , &
               dv   => block%dtend(1)%dv, &
               u    => dstate%u_lon     , &
               v    => dstate%v_lat     )
    select case (vor_damp_order)
    case (2)
      call fill_halo(vor, east_halo=.false., north_halo=.false.)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u%d(i,j,k) = u%d(i,j,k) - dt * cx(j,k) * (vor%d(i,j,k) - vor%d(i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dv%d(i,j,k) = dt * cy(j,k) * (vor%d(i,j,k) - vor%d(i-1,j,k)) / mesh%le_lat(j)
          end do
        end do
      end do
      call fill_halo(dv, south_halo=.false., north_halo=.false.)
      call filter_on_lat_edge(block%small_filter, dv%d)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = v%d(i,j,k) + dv%d(i,j,k)
          end do
        end do
      end do
    end select
    call fill_halo(u)
    call fill_halo(v)
    end associate

  end subroutine vor_damp_run

end module vor_damp_mod
