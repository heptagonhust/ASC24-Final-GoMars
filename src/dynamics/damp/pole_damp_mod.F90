module pole_damp_mod

  use namelist_mod
  use time_mod
  use math_mod
  use block_mod
  use filter_mod
  use laplace_damp_mod
  use operators_mod
  use parallel_mod

  implicit none

  private

  public pole_damp_run

contains

  subroutine pole_damp_run(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) c, tmp(global_mesh%full_nlev)

    associate (mesh    => block%mesh  , &
               dmg     => dstate%dmg  , &
               pt      => dstate%pt   , &
               qv      => dstate%qv   , &
               u_lon   => dstate%u_lon, &
               v_lat   => dstate%v_lat)
    if (use_pole_damp .and. baroclinic) then
      ! c = 1.0e12_r8
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt(i,j,k) = pt(i,j,k) * dmg(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      ! call laplace_damp_on_cell(block%filter_mesh, block%filter_halo, 4, pt, lon_coef=decay_from_pole, coef=c)
      call filter_on_cell(block%small_filter, pt)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            pt(i,j,k) = pt(i,j,k) / dmg(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
      ! ----------------------------------------------------------------------
      if (time_is_alerted('moist')) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qv(i,j,k) = qv(i,j,k) * dmg(i,j,k)
            end do
          end do
        end do
        call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true.)
        ! call laplace_damp_on_cell(block%filter_mesh, block%filter_halo, 4, qv, lon_coef=decay_from_pole, coef=c)
        call filter_on_cell(block%small_filter, qv)
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              qv(i,j,k) = qv(i,j,k) / dmg(i,j,k)
            end do
          end do
        end do
        call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
      end if
      ! ----------------------------------------------------------------------
      ! call laplace_damp_on_lon_edge(block%filter_mesh, block%filter_halo, 4, u_lon, lon_coef=decay_from_pole, coef=c)
      ! call laplace_damp_on_lat_edge(block%filter_mesh, block%filter_halo, 4, v_lat, lon_coef=decay_from_pole, coef=c)
      ! ----------------------------------------------------------------------
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        c = exp_two_values(0.05_r8, 0.0_r8, 90.0_r8, 85.0_r8, real(abs(mesh%full_lat_deg(j)), r8))
        do i = mesh%half_ids, mesh%half_ide
          tmp = u_lon(i,j,1:mesh%full_nlev)
          do k = mesh%full_kds + 1, mesh%full_kde - 1
            u_lon(i,j,k) = c * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c) * tmp(k)
          end do
        end do
      end do
      call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      do j = mesh%half_jds, mesh%half_jde
        c = exp_two_values(0.05_r8, 0.0_r8, 90.0_r8, 85.0_r8, real(abs(mesh%half_lat_deg(j)), r8))
        do i = mesh%full_ids, mesh%full_ide
          tmp = v_lat(i,j,1:mesh%full_nlev)
          do k = mesh%full_kds + 1, mesh%full_kde - 1
            v_lat(i,j,k) = c * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c) * tmp(k)
          end do
        end do
      end do
      call fill_halo(block%halo, v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
      ! ----------------------------------------------------------------------
    end if

    if (nudge_pole_v) then
      c = nudge_pole_v_coef
      do j = mesh%half_jms, mesh%half_jme
        if (mesh%is_south_pole(j)) then
          v_lat(:,j,:) = (1 - c) * v_lat(:,j,:) + c * v_lat(:,j+1,:)
        else if (mesh%is_north_pole(j+1)) then
          v_lat(:,j,:) = (1 - c) * v_lat(:,j,:) + c * v_lat(:,j-1,:)
        end if
      end do
    end if
    end associate

  end subroutine pole_damp_run

end module pole_damp_mod
