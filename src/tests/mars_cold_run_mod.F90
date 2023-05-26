module mars_cold_run_mod

  use flogger
  use namelist_mod, only: topo_file, use_topo_smooth, ptop
  use const_mod
  use latlon_parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use topo_mod
  use operators_mod

  implicit none

  private

  public mars_cold_run_set_ic

  real(r8), parameter :: t0   = 170.0_r8   ! K
  real(r8), parameter :: ps0  = 610.0_r8   ! Pa

contains

  subroutine mars_cold_run_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) min_lon, max_lon, min_lat, max_lat
    integer i, j, k

    associate (mesh   => block%mesh            , &
               u      => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v_lat , &
               t      => block%dstate(1)%t     , &
               pt     => block%dstate(1)%pt    , &
               mg     => block%dstate(1)%mg    , &
               mgs    => block%dstate(1)%mgs   , &
               phs    => block%dstate(1)%phs   , &
               gzs    => block%static%gzs)
    min_lon = mesh%full_lon_deg(mesh%full_ids-1)
    max_lon = mesh%full_lon_deg(mesh%full_ide+1)
    min_lat = mesh%full_lat_deg(max(1, mesh%full_jds-1))
    max_lat = mesh%full_lat_deg(min(global_mesh%full_nlat, mesh%full_jde+1))
    call topo_read(min_lon, max_lon, min_lat, max_lat)
    call topo_regrid(block)
    if (use_topo_smooth) then
      call topo_smooth(block)
    end if

    u = 0
    v = 0
    t = t0
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        mgs(i,j) = ps0 * exp(-gzs(i,j) / (Rd * t0)) - ptop
      end do
    end do
    call fill_halo(block%halo, mgs, full_lon=.true., full_lat=.true.)

    call calc_mg(block, block%dstate(1))

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          pt(i,j,k) = modified_potential_temperature(t(i,j,k), mg(i,j,k), 0.0_r8)
        end do
      end do
    end do
    call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    phs = mgs
    end associate

  end subroutine mars_cold_run_set_ic

end module mars_cold_run_mod
