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

  subroutine pole_damp_run(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    if (baroclinic) then
      state%pt = state%pt * state%m
      call filter_on_cell(block%small_filter_phs, state%phs)
      call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)
      call calc_ph(block, state)
      call calc_m (block, state)
      call filter_on_cell(block%small_filter_pt, state%pt)
      state%pt = state%pt / state%m
      call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      call filter_on_lon_edge(block%small_filter_uv, state%u_lon)
      call fill_halo(block, state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      call filter_on_lat_edge(block%small_filter_uv, state%v_lat)
      call fill_halo(block, state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
    else
      ! call filter_on_cell(block%small_filter_phs, state%gz)
      ! call fill_halo(block, state%gz, full_lon=.true., full_lat=.true.)
      ! call calc_m (block, state)
    end if

  end subroutine pole_damp_run

end module pole_damp_mod
