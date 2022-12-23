module pole_damp_mod

  use namelist_mod
  use block_mod
  use filter_mod
  use operators_mod
  use parallel_mod

  implicit none

  private

  public pole_damp_run

contains

  subroutine pole_damp_run(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    if (baroclinic) then
      dstate%pt = dstate%pt * dstate%m
      call filter_on_cell(block%small_filter, dstate%phs)
      call fill_halo(block%halo, dstate%phs, full_lon=.true., full_lat=.true.)
      call calc_ph(block, dstate)
      call calc_m (block, dstate)
      call filter_on_cell(block%small_filter, dstate%pt)
      dstate%pt = dstate%pt / dstate%m
      call fill_halo(block%halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      call filter_on_lon_edge(block%small_filter, dstate%u_lon)
      call fill_halo(block%halo, dstate%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      call filter_on_lat_edge(block%small_filter, dstate%v_lat)
      call fill_halo(block%halo, dstate%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
    else
      ! call filter_on_cell(block%small_filter, dstate%gz)
      ! call fill_halo(block%halo, dstate%gz, full_lon=.true., full_lat=.true.)
      ! call calc_m (block, dstate)
    end if

  end subroutine pole_damp_run

end module pole_damp_mod
