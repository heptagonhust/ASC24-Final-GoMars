module div_damp_mod

  use flogger
  use string
  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use process_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: c2_lon(:,:), c4_lon(:,:)
  real(r8), allocatable :: c2_lat(:,:), c4_lat(:,:)

contains

  subroutine div_damp_init()

    real(r8) r, c0
    integer j, k

    if (.not. use_div_damp) return

    call div_damp_final()

    allocate(c2_lon(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(c2_lat(global_mesh%half_nlat,global_mesh%full_nlev))
    allocate(c4_lon(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(c4_lat(global_mesh%half_nlat,global_mesh%full_nlev))

    c0 = merge(0.0_r8, 1.0_r8, div_damp_order > 2)
    r = 1
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
        c2_lon(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**(r - 1) * &
          exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
          exp_two_values(div_damp_pole, c0, 90.0_r8, 80.0_r8, abs(global_mesh%full_lat_deg(j))) * &
          global_mesh%le_lon(j) * global_mesh%de_lon(j)
      end do
    end do
    r = 1
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%half_jds, global_mesh%half_jde
        c2_lat(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**(r - 1) * &
          exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
          exp_two_values(div_damp_pole, c0, 90.0_r8, 80.0_r8, abs(global_mesh%half_lat_deg(j))) * &
          global_mesh%le_lat(j) * global_mesh%de_lat(j)
      end do
    end do
    r = 2
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
        c4_lon(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**(r - 2) * &
          global_mesh%le_lon(j)**2 * global_mesh%de_lon(j)**2
      end do
    end do
    r = 2
    do k = global_mesh%full_kds, global_mesh%full_kde
      do j = global_mesh%half_jds, global_mesh%half_jde
        c4_lat(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**(r - 2) * &
          global_mesh%le_lat(j)**2 * global_mesh%de_lat(j)**2
      end do
    end do

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c2_lon)) deallocate(c2_lon)
    if (allocated(c2_lat)) deallocate(c2_lat)
    if (allocated(c4_lon)) deallocate(c4_lon)
    if (allocated(c4_lat)) deallocate(c4_lat)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(8), intent(in) :: dt

    integer i, j, k

    call calc_div(block, dstate)

    associate (mesh => block%mesh         , &
               div  => block%aux%div      , &
               div2 => block%aux%div2     , &
               dudt => block%aux%dudt_damp, &
               dvdt => block%aux%dvdt_damp, &
               u    => dstate%u_lon       , &
               v    => dstate%v_lat       )
    if (div_damp_order >= 2) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u(i,j,k) = u(i,j,k) + c2_lon(j,k) * (div(i+1,j,k) - div(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v(i,j,k) = v(i,j,k) + c2_lat(j,k) * (div(i,j+1,k) - div(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    end if
    if (div_damp_order == 4) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            u(i,j,k) = u(i,j,k) - c4_lon(j,k) * (div2(i+1,j,k) - div2(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            v(i,j,k) = v(i,j,k) - c4_lat(j,k) * (div2(i,j+1,k) - div2(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    end if
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine div_damp_run

end module div_damp_mod
