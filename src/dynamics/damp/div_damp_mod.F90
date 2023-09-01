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

  real(r8), allocatable :: c_lon(:,:)
  real(r8), allocatable :: c_lat(:,:)

contains

  subroutine div_damp_init()

    integer j, k

    if (.not. use_div_damp) return

    call div_damp_final()

    allocate(c_lon(global_mesh%full_nlat,global_mesh%full_nlev))
    allocate(c_lat(global_mesh%half_nlat,global_mesh%full_nlev))

    select case (div_damp_order)
    case (2)
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
          if (baroclinic) then
            c_lon(j,k) = div_damp_coef2 * &
              exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
              global_mesh%le_lon(j) * global_mesh%de_lon(j)
          else
            c_lon(j,k) = div_damp_coef2 * global_mesh%le_lon(j) * global_mesh%de_lon(j)
          end if
        end do
      end do
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%half_jds, global_mesh%half_jde
          if (baroclinic) then
            c_lat(j,k) = div_damp_coef2 * &
              exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
              global_mesh%le_lat(j) * global_mesh%de_lat(j)
          else
            c_lat(j,k) = div_damp_coef2 * global_mesh%le_lat(j) * global_mesh%de_lat(j)
          end if
        end do
      end do
    case (4)
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
          c_lon(j,k) = div_damp_coef4 * &
            ! exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
            global_mesh%le_lon(j)**2 * global_mesh%de_lon(j)**2
        end do
      end do
      do k = global_mesh%full_kds, global_mesh%full_kde
        do j = global_mesh%half_jds, global_mesh%half_jde
          c_lat(j,k) = div_damp_coef4 * &
            ! exp_two_values(div_damp_top, 1.0_r8, 1.0_r8, real(div_damp_k0, r8), real(k, r8)) * &
            global_mesh%le_lat(j)**2 * global_mesh%de_lat(j)**2
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_str(div_damp_order)) // '!')
    end select

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    call calc_div(block, dstate)

    associate (mesh => block%mesh        , &
               div  => block%aux%div     , &
               div2 => block%aux%div2    , &
               dudt => block%aux%dudt_div, &
               dvdt => block%aux%dvdt_div, &
               u    => dstate%u_lon      , &
               v    => dstate%v_lat      )
    select case (div_damp_order)
    case (2)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            dudt(i,j,k) = c_lon(j,k) * (div(i+1,j,k) - div(i,j,k)) / mesh%de_lon(j) / dt
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dvdt(i,j,k) = c_lat(j,k) * (div(i,j+1,k) - div(i,j,k)) / mesh%de_lat(j) / dt
          end do
        end do
      end do
    case (4)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            dudt(i,j,k) = -c_lon(j,k) * (div2(i+1,j,k) - div2(i,j,k)) / mesh%de_lon(j) / dt
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dvdt(i,j,k) = -c_lat(j,k) * (div2(i,j+1,k) - div2(i,j,k)) / mesh%de_lat(j) / dt
          end do
        end do
      end do
    end select
    end associate

  end subroutine div_damp_run

end module div_damp_mod
