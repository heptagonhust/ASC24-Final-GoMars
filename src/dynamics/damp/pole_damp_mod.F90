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

    allocate(c_lon(global_mesh%full_jms:global_mesh%full_jme)); c_lon = 0
    allocate(c_lat(global_mesh%half_jms:global_mesh%half_jme)); c_lat = 0

    do j = global_mesh%full_jms, global_mesh%full_jme
      if (j > global_mesh%full_jds .and. j < global_mesh%full_jde) then
        c_lon(j) = exp_two_values(1.0_r8, 0.0_r8               , &
          real(abs(global_mesh%full_lat_deg(2)), r8)           , &
          pole_damp_lat0                                       , &
          real(abs(global_mesh%full_lat_deg(j)), r8))
        c_lon(j) = merge(0.0_r8, c_lon(j), c_lon(j) < 1.0e-15)
      end if
    end do

    do j = global_mesh%half_jms, global_mesh%half_jme
      if (j >= global_mesh%half_jds .and. j <= global_mesh%half_jde) then
        c_lat(j) = exp_two_values(1.0_r8, 0.0_r8               , &
          real(abs(global_mesh%half_lat_deg(1)), r8)           , &
          pole_damp_lat0                                       , &
          real(abs(global_mesh%half_lat_deg(j)), r8))
        c_lat(j) = merge(0.0_r8, c_lat(j), c_lat(j) < 1.0e-15)
      end if
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

    real(r8) tmp(block%mesh%full_ims:block%mesh%full_ime, &
                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8) wgt1, wgt2, c
    integer i, j, k

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, &
               v    => dstate%v_lat, &
               dmg  => dstate%dmg  , &
               pt   => dstate%pt   )
    do j = mesh%half_jms, mesh%half_jme
      if (c_lat(j) > 0) then
        tmp = v(:,j,:)
        do i = mesh%full_ims + 1, mesh%full_ime - 1
          v(i,j,:) = c_lat(j) * (tmp(i,:) + 0.25_r8 * (tmp(i-1,:) + tmp(i+1,:) - 2 * tmp(i,:))) + &
                     (1 - c_lat(j)) * tmp(i,:)
        end do
      end if
    end do
    ! This nudging of polar v helps to keep the flow neat around the poles.
    ! NOTE: DO NOT REMOVE IT!
    c = 0.2_r8
    do j = mesh%half_jms, mesh%half_jme
      if (mesh%is_south_pole(j)) then
        v(:,j,:) = (1 - c) * v(:,j,:) + c * v(:,j+1,:)
      else if (mesh%is_north_pole(j+1)) then
        v(:,j,:) = (1 - c) * v(:,j,:) + c * v(:,j-1,:)
      end if
    end do
    ! do j = mesh%full_jms, mesh%full_jme
    !   if (c_lon(j) > 0) then
    !     tmp = u(:,j,:)
    !     do k = mesh%full_kds + 1, mesh%full_kde - 1
    !       wgt1 = mesh%half_dlev(k+1) / (mesh%half_dlev(k) + mesh%half_dlev(k+1))
    !       wgt2 = 1 - wgt1
    !       u(:,j,k) = tmp(:,k) - c_lon(j) * (tmp(:,k) - (wgt1 * tmp(:,k-1) + wgt2 * tmp(:,k+1)))
    !     end do
    !     k = mesh%full_kds
    !     wgt1 = (mesh%half_dlev(k+1) + mesh%half_dlev(k+2)) / mesh%half_dlev(k+2)
    !     wgt2 = -mesh%half_dlev(k+1) / mesh%half_dlev(k+2)
    !     u(:,j,k) = tmp(:,k) - c_lon(j) * (tmp(:,k) - (wgt1 * tmp(:,k+1) + wgt2 * tmp(:,k+2)))
    !     k = mesh%full_kde
    !     wgt1 = (mesh%half_dlev(k) + mesh%half_dlev(k-1)) / mesh%half_dlev(k-1)
    !     wgt2 = -mesh%half_dlev(k) / mesh%half_dlev(k-1)
    !     u(:,j,k) = tmp(:,k) - c_lon(j) * (tmp(:,k) - (wgt1 * tmp(:,k-1) + wgt2 * tmp(:,k-2)))
    !   end if
    ! end do
    ! do j = mesh%half_jds, mesh%half_jde
    !   if (c_lat(j) > 0) then
    !     tmp = v(:,j,:)
    !     do k = mesh%full_kds + 1, mesh%full_kde - 1
    !       wgt1 = mesh%half_dlev(k+1) / (mesh%half_dlev(k) + mesh%half_dlev(k+1))
    !       wgt2 = 1 - wgt1
    !       v(:,j,k) = tmp(:,k) - c_lat(j) * (tmp(:,k) - (wgt1 * tmp(:,k-1) + wgt2 * tmp(:,k+1)))
    !     end do
    !     k = mesh%full_kds
    !     wgt1 = (mesh%half_dlev(k+1) + mesh%half_dlev(k+2)) / mesh%half_dlev(k+2)
    !     wgt2 = -mesh%half_dlev(k+1) / mesh%half_dlev(k+2)
    !     v(:,j,k) = tmp(:,k) - c_lat(j) * (tmp(:,k) - (wgt1 * tmp(:,k+1) + wgt2 * tmp(:,k+2)))
    !     k = mesh%full_kde
    !     wgt1 = (mesh%half_dlev(k) + mesh%half_dlev(k-1)) / mesh%half_dlev(k-1)
    !     wgt2 = -mesh%half_dlev(k) / mesh%half_dlev(k-1)
    !     v(:,j,k) = tmp(:,k) - c_lat(j) * (tmp(:,k) - (wgt1 * tmp(:,k-1) + wgt2 * tmp(:,k-2)))
    !   end if
    ! end do
    end associate

  end subroutine pole_damp_run

end module pole_damp_mod
