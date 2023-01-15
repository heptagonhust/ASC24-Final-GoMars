module time_schemes_mod

  use flogger
  use const_mod
  use math_mod
  use namelist_mod
  use block_mod
  use operators_mod
  use parallel_mod
  use filter_mod

  implicit none

  private

  public time_scheme_init
  public time_scheme_final
  public time_integrator
  public update_state
  public space_operators_interface

  interface
    subroutine space_operators_interface(block, old_state, star_state, new_state, tend1, tend2, dt, pass)
      import block_type, dstate_type, dtend_type, r8
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_state
      type(dstate_type), intent(inout) :: star_state
      type(dstate_type), intent(inout) :: new_state
      type(dtend_type ), intent(inout) :: tend1
      type(dtend_type ), intent(in   ) :: tend2
      real(r8), intent(in) :: dt
      integer, intent(in) :: pass
    end subroutine space_operators_interface

    subroutine step_interface(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)
      import space_operators_interface, block_type, dstate_type, dtend_type, r8
      procedure(space_operators_interface) space_operators
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_state
      type(dstate_type), intent(inout) :: star_state
      type(dstate_type), intent(inout) :: new_state
      type(dtend_type ), intent(inout) :: tend1
      type(dtend_type ), intent(inout) :: tend2
      real(r8), intent(in) :: dt
    end subroutine step_interface

    subroutine time_integrator_interface(space_operators, block, old, new, dt)
      import block_type, dtend_type, dstate_type, space_operators_interface, r8
      procedure(space_operators_interface) space_operators
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      real(r8), intent(in) :: dt
    end subroutine time_integrator_interface
  end interface

  procedure(step_interface), pointer :: step
  procedure(time_integrator_interface), pointer :: time_integrator

contains

  subroutine time_scheme_init()

    call time_scheme_final()

    select case (time_scheme)
    case ('euler')
      time_integrator => euler
    case ('pc2')
      time_integrator => pc2
    case ('wrfrk3')
      time_integrator => wrfrk3
    case default
      time_integrator => pc2
    end select

    step => step_forward_backward

  end subroutine time_scheme_init

  subroutine time_scheme_final()

  end subroutine time_scheme_final

  subroutine step_all(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type ), intent(inout) :: tend1
    type(dtend_type ), intent(inout) :: tend2
    real(r8), intent(in) :: dt

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, all_pass)
    call update_state(block, tend1, old_state, new_state, dt)

  end subroutine step_all

  subroutine step_forward_backward(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type ), intent(inout) :: tend1
    type(dtend_type ), intent(inout) :: tend2
    real(r8), intent(in) :: dt

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, forward_pass)
    call update_state(block, tend1, old_state, new_state, dt)
    call space_operators(block, old_state, star_state, new_state, tend2, tend1, dt, backward_pass)
    call update_state(block, tend2, old_state, new_state, dt)

  end subroutine step_forward_backward

  subroutine update_state(block, dtend, old_state, new_state, dt)

    type(block_type ), intent(inout) :: block
    type(dtend_type ), intent(inout) :: dtend
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: new_state
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) c, tmp(global_mesh%full_nlev)

    associate (mesh => block%mesh)
    if (baroclinic) then
      if (dtend%update_phs) then
        ! ----------------------------------------------------------------------
        call fill_halo(block%halo, dtend%dphs, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, dtend%dphs)
        ! ----------------------------------------------------------------------
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            new_state%phs(i,j) = old_state%phs(i,j) + dt * dtend%dphs(i,j)
          end do
        end do
        call fill_halo(block%halo, new_state%phs, full_lon=.true., full_lat=.true.)
        call calc_ph(block, new_state)
        call calc_m (block, new_state)
      else if (dtend%copy_phs) then
        new_state%phs    = old_state%phs
        new_state%ph_lev = old_state%ph_lev
        new_state%ph     = old_state%ph
        new_state%m      = old_state%m
      end if

      if (dtend%update_pt) then
        if (.not. dtend%update_phs .and. .not. dtend%copy_phs .and. proc%is_root()) call log_error('Mass is not updated or copied!')
        ! ----------------------------------------------------------------------
        call fill_halo(block%halo, dtend%dpt, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, dtend%dpt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_state%pt(i,j,k) = (old_state%pt(i,j,k) * old_state%m(i,j,k) + dt * dtend%dpt(i,j,k)) / new_state%m(i,j,k)
            end do
          end do
        end do
        call fill_halo(block%halo, new_state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      else if (dtend%copy_pt) then
        new_state%pt = old_state%pt
      end if
    else
      if (dtend%update_gz) then
        ! ----------------------------------------------------------------------
        call fill_halo(block%halo, dtend%dgz, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, dtend%dgz)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_state%gz(i,j,k) = old_state%gz(i,j,k) + dt * dtend%dgz(i,j,k)
            end do
          end do
        end do
        call fill_halo(block%halo, new_state%gz, full_lon=.true., full_lat=.true.)
        call calc_m(block, new_state)
      else if (dtend%copy_gz) then
        new_state%gz = old_state%gz
      end if
    end if

    if (dtend%update_u .and. dtend%update_v) then
      ! ----------------------------------------------------------------------
      call fill_halo(block%halo, dtend%du, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
      call filter_on_lon_edge(block%big_filter, dtend%du)
      call fill_halo(block%halo, dtend%dv, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
      call filter_on_lat_edge(block%big_filter, dtend%dv)
      ! ----------------------------------------------------------------------
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            new_state%u_lon(i,j,k) = old_state%u_lon(i,j,k) + dt * dtend%du(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            new_state%v_lat(i,j,k) = old_state%v_lat(i,j,k) + dt * dtend%dv(i,j,k)
          end do
        end do
      end do
      ! ----------------------------------------------------------------------
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        c = exp_two_values(0.1_r8, 0.0_r8, 90.0_r8, 85.0_r8, abs(mesh%full_lat_deg(j)))
        do i = mesh%half_ids, mesh%half_ide
          tmp = new_state%u_lon(i,j,1:mesh%full_nlev)
          do k = mesh%full_kds + 1, mesh%full_kde - 1
            new_state%u_lon(i,j,k) = c * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c) * tmp(k)
          end do
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        c = exp_two_values(0.1_r8, 0.0_r8, 90.0_r8, 85.0_r8, abs(mesh%half_lat_deg(j)))
        do i = mesh%full_ids, mesh%full_ide
          tmp = new_state%v_lat(i,j,1:mesh%full_nlev)
          do k = mesh%full_kds + 1, mesh%full_kde - 1
            new_state%v_lat(i,j,k) = c * (tmp(k-1) + tmp(k+1)) + (1 - 2 * c) * tmp(k)
          end do
        end do
      end do
      ! ----------------------------------------------------------------------
      call fill_halo(block%halo, new_state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      call fill_halo(block%halo, new_state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
    end if
    end associate

  end subroutine update_state

  subroutine pc2(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt / 2.0_r8)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend(old), dtend(new), dt / 2.0_r8)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend(old), dtend(new), dt         )
    end associate

  end subroutine pc2

  subroutine wrfrk3(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt / 3.0_r8)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend(old), dtend(new), dt / 2.0_r8)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend(old), dtend(new), dt         )
    end associate

  end subroutine wrfrk3

  subroutine euler(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt)
    end associate

  end subroutine euler

end module time_schemes_mod
