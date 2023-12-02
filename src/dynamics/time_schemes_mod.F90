! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module implements time integration schemes, including PC2 and WRFRK3,
!   and with forward-backward technique.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
!   - Jianghao Li
! ==============================================================================

module time_schemes_mod

  use flogger
  use const_mod
  use namelist_mod
  use block_mod
  use operators_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use filter_mod

  implicit none

  private

  public time_scheme_init
  public time_scheme_final
  public time_integrator
  public update_state
  public space_operators_interface

  interface
    subroutine space_operators_interface(block, old_state, star_state, new_state, tend1, tend2, dt, pass, substep)
      import block_type, dstate_type, dtend_type, r8
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_state
      type(dstate_type), intent(inout) :: star_state
      type(dstate_type), intent(inout) :: new_state
      type(dtend_type ), intent(inout) :: tend1
      type(dtend_type ), intent(in   ) :: tend2
      real(r8), intent(in) :: dt
      integer, intent(in) :: pass
      integer, intent(in) :: substep
    end subroutine space_operators_interface

    subroutine step_interface(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt, substep)
      import space_operators_interface, block_type, dstate_type, dtend_type, r8
      procedure(space_operators_interface) space_operators
      type(block_type ), intent(inout) :: block
      type(dstate_type), intent(in   ) :: old_state
      type(dstate_type), intent(inout) :: star_state
      type(dstate_type), intent(inout) :: new_state
      type(dtend_type ), intent(inout) :: tend1
      type(dtend_type ), intent(inout) :: tend2
      real(r8), intent(in) :: dt
      integer, intent(in) :: substep
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
      total_substeps = 1
    case ('rk2')
      time_integrator => rk2
      total_substeps = 2
    case ('pc2')
      time_integrator => pc2
      total_substeps = 3
    case ('wrfrk3')
      time_integrator => wrfrk3
      total_substeps = 3
    case default
      time_integrator => pc2
      total_substeps = 3
    end select

    step => step_forward_backward

  end subroutine time_scheme_init

  subroutine time_scheme_final()

  end subroutine time_scheme_final

  subroutine step_all(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt, substep)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type ), intent(inout) :: tend1
    type(dtend_type ), intent(inout) :: tend2
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, all_pass, substep)
    call update_state(block, tend1, old_state, new_state, dt)

  end subroutine step_all

  subroutine step_forward_backward(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt, substep)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type ), intent(inout) :: tend1
    type(dtend_type ), intent(inout) :: tend2
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, forward_pass, substep)
    call update_state(block, tend1, old_state, new_state, dt)
    call space_operators(block, old_state, star_state, new_state, tend2, tend1, dt, backward_pass, substep)
    call update_state(block, tend2, old_state, new_state, dt)

  end subroutine step_forward_backward

  subroutine update_state(block, dtend, old_state, new_state, dt)

    type(block_type ), intent(inout) :: block
    type(dtend_type ), intent(inout) :: dtend
    type(dstate_type), intent(in   ) :: old_state
    type(dstate_type), intent(inout) :: new_state
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh       => block%mesh, &
               dmgsdt     => dtend%dmgs, &
               dgzdt      => dtend%dgz , &
               dptdt      => dtend%dpt , &
               dudt       => dtend%du  , &
               dvdt       => dtend%dv  )
    if (baroclinic) then
      if (dtend%update_mgs) then
        ! ----------------------------------------------------------------------
        call fill_halo(dmgsdt, south_halo=.false., north_halo=.false.)
        call filter_run(block%big_filter, dmgsdt)
        ! ----------------------------------------------------------------------
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            new_state%mgs%d(i,j) = old_state%mgs%d(i,j) + dt * dmgsdt%d(i,j)
          end do
        end do
        call fill_halo(new_state%mgs)
        call calc_mg (block, new_state)
        call calc_dmg(block, new_state)
        call calc_ph (block, new_state)
      else if (dtend%copy_mgs) then ! FIXME: Do we still need copy?
        new_state%mgs   %d = old_state%mgs   %d
        new_state%mg_lev%d = old_state%mg_lev%d
        new_state%mg    %d = old_state%mg    %d
        new_state%dmg   %d = old_state%dmg   %d
        new_state%ph_lev%d = old_state%ph_lev%d
        new_state%ph    %d = old_state%ph    %d
      end if

      if (dtend%update_pt) then
        if (.not. dtend%update_mgs .and. .not. dtend%copy_mgs .and. proc%is_root()) call log_error('Mass is not updated or copied!')
        ! ----------------------------------------------------------------------
        call fill_halo(dptdt, south_halo=.false., north_halo=.false.)
        call filter_run(block%big_filter, dptdt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_state%pt%d(i,j,k) = (old_state%pt%d(i,j,k) * old_state%dmg%d(i,j,k) + dt * dptdt%d(i,j,k)) / new_state%dmg%d(i,j,k)
            end do
          end do
        end do
        call fill_halo(new_state%pt)
      else if (dtend%copy_pt) then
        new_state%pt%d = old_state%pt%d
      end if
    else
      if (dtend%update_gz) then
        ! ----------------------------------------------------------------------
        call fill_halo(dgzdt, south_halo=.false., north_halo=.false.)
        call filter_run(block%big_filter, dgzdt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              new_state%gz%d(i,j,k) = old_state%gz%d(i,j,k) + dt * dgzdt%d(i,j,k)
            end do
          end do
        end do
        call fill_halo(new_state%gz)
        call calc_dmg(block, new_state)
      else if (dtend%copy_gz) then
        new_state%gz = old_state%gz
      end if
    end if

    if (dtend%update_u .and. dtend%update_v) then
      ! ------------------------------------------------------------------------
      call fill_halo(dudt, south_halo=.false., north_halo=.false.)
      call filter_run(block%big_filter, dudt)
      call fill_halo(dvdt, south_halo=.false., north_halo=.false.)
      call filter_run(block%big_filter, dvdt)
      ! ------------------------------------------------------------------------
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            new_state%u_lon%d(i,j,k) = old_state%u_lon%d(i,j,k) + dt * dudt%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            new_state%v_lat%d(i,j,k) = old_state%v_lat%d(i,j,k) + dt * dvdt%d(i,j,k)
          end do
        end do
      end do
      call fill_halo(new_state%u_lon)
      call fill_halo(new_state%v_lat)
    end if
    end associate

  end subroutine update_state

  subroutine rk2(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type ), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt / 2.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(new), dtend(old), dtend(new), dt         , 2)
    end associate

  end subroutine rk2

  subroutine pc2(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt / 2.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend(old), dtend(new), dt / 2.0_r8, 2)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend(old), dtend(new), dt         , 3)
    end associate

  end subroutine pc2

  subroutine wrfrk3(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt / 3.0_r8, 1)
    call step(space_operators, block, dstate(old), dstate(new), dstate(3  ), dtend(old), dtend(new), dt / 2.0_r8, 2)
    call step(space_operators, block, dstate(old), dstate(3  ), dstate(new), dtend(old), dtend(new), dt         , 3)
    end associate

  end subroutine wrfrk3

  subroutine euler(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(r8), intent(in) :: dt

    associate (dstate => block%dstate, dtend => block%dtend)
    call step(space_operators, block, dstate(old), dstate(old), dstate(new), dtend(old), dtend(new), dt, 1)
    end associate

  end subroutine euler

end module time_schemes_mod
