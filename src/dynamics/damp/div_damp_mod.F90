! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module div_damp_mod

  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use operators_mod
  use filter_mod

  implicit none

  private

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: cx(:,:), cy(:,:), cx_pole(:,:)

contains

  subroutine div_damp_init()

    real(r8) r, lat0
    integer j, k

    if (.not. use_div_damp) return

    call div_damp_final()

    allocate(cx(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(cy(global_mesh%half_nlat,global_mesh%full_nlev))
    allocate(cx_pole(global_mesh%full_nlat,global_mesh%full_nlev))

    select case (div_damp_order)
    case (2)
      r = 1
      lat0 = abs(global_mesh%full_lat_deg(2))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
          cx(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**(r - 1) * &
            exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
            global_mesh%le_lon(j) * global_mesh%de_lon(j) / dt_dyn
          cx_pole(j,k) = cx(j,k) * &
            exp_two_values(div_damp_pole, 0.0_r8, lat0, div_damp_lat0, abs(global_mesh%full_lat_deg(j)))
        end do
      end do
      lat0 = abs(global_mesh%half_lat_deg(1))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%half_jds, global_mesh%half_jde
          cy(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**(r - 1) * &
            exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
            global_mesh%le_lat(j) * global_mesh%de_lat(j) / dt_dyn
        end do
      end do
    case (4)
      r = 2
      lat0 = abs(global_mesh%full_lat_deg(2))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
          cx(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**(r - 2) * &
            global_mesh%le_lon(j)**2 * global_mesh%de_lon(j)**2
        end do
      end do
      lat0 = abs(global_mesh%half_lat_deg(1))
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%half_jds, global_mesh%half_jde
          cy(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**(r - 2) * &
            global_mesh%le_lat(j)**2 * global_mesh%de_lat(j)**2
        end do
      end do
    end select

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(cx)) deallocate(cx)
    if (allocated(cy)) deallocate(cy)
    if (allocated(cx_pole)) deallocate(cx_pole)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(8), intent(in) :: dt

    integer i, j, k

    call calc_div(block, dstate)

    associate (mesh => block%mesh       , &
               div  => block%aux%div    , &
               div2 => block%aux%div2   , &
               du   => block%dtend(1)%du, &
               u    => dstate%u_lon     , &
               v    => dstate%v_lat     )
    select case (div_damp_order)
    case (2)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u%d(i,j,k) = u%d(i,j,k) + dt * cx(j,k) * (div%d(i+1,j,k) - div%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = v%d(i,j,k) + dt * cy(j,k) * (div%d(i,j+1,k) - div%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      ! ------------------------------------------------------------------------
      if (div_damp_pole > 0) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              du%d(i,j,k) = dt * cx_pole(j,k) * (div%d(i+1,j,k) - div%d(i,j,k)) / mesh%de_lon(j)
            end do
          end do
        end do
        call fill_halo(du, south_halo=.false., north_halo=.false.)
        call filter_run(block%small_filter, du)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              u%d(i,j,k) = u%d(i,j,k) + du%d(i,j,k)
            end do
          end do
        end do
      end if
      ! ------------------------------------------------------------------------
    case (4)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            du%d(i,j,k) = -cx(j,k) * (div2%d(i+1,j,k) - div2%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      call fill_halo(du, south_halo=.false., north_halo=.false.)
      call filter_run(block%small_filter, du)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u%d(i,j,k) = u%d(i,j,k) + du%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = v%d(i,j,k) - cy(j,k) * (div2%d(i,j+1,k) - div2%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    end select
    call fill_halo(u)
    call fill_halo(v)
    end associate

  end subroutine div_damp_run

end module div_damp_mod
