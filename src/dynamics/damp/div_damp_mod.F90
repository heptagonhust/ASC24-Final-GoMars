module div_damp_mod

  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: c2(:,:), c4(:,:)

contains

  subroutine div_damp_init()

    real(r8) r, c0
    integer j, k

    if (.not. use_div_damp) return

    call div_damp_final()

    allocate(c2(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(c4(global_mesh%half_nlat,global_mesh%full_nlev))

    c0 = merge(0.0_r8, 1.0_r8, div_damp_order > 2)
    r = 1
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
        c2(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**(r - 1) * &
          exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
          exp_two_values(div_damp_pole, c0, 90.0_r8, 80.0_r8, abs(global_mesh%full_lat_deg(j))) * &
          global_mesh%le_lon(j) * global_mesh%de_lon(j)
      end do
    end do
    r = 2
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
        c4(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**(r - 2) * &
          global_mesh%le_lon(j)**2 * global_mesh%de_lon(j)**2
      end do
    end do

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c2)) deallocate(c2)
    if (allocated(c4)) deallocate(c4)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(8), intent(in) :: dt

    integer i, j, k

    call calc_div(block, dstate)

    associate (mesh => block%mesh    , &
               div  => block%aux%div , &
               div2 => block%aux%div2, &
               u    => dstate%u_lon  , &
               v    => dstate%v_lat  )
    if (div_damp_order >= 2) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u(i,j,k) = u(i,j,k) + c2(j,k) * (div(i+1,j,k) - div(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v(i,j,k) = v(i,j,k) + (c2(j+1,k) * div(i,j+1,k) - c2(j,k) * div(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    end if
    if (div_damp_order == 4) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u(i,j,k) = u(i,j,k) - c4(j,k) * (div2(i+1,j,k) - div2(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v(i,j,k) = v(i,j,k) - (c4(j+1,k) * div2(i,j+1,k) - c4(j,k) * div2(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    end if
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine div_damp_run

end module div_damp_mod
