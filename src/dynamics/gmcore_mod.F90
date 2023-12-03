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
!   This is the main module of GMCORE, which provides initialization, run,
!   finalization subroutines.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
!   - Jianghao Li
! ==============================================================================

module gmcore_mod

  use mpi
  use flogger
  use string
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use process_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use restart_mod
  use block_mod
  use vert_coord_mod
  use time_schemes_mod
  use operators_mod
  use tracer_mod
  use interp_mod
  use gas_mod
  use adv_mod
  use pgf_mod
  use damp_mod
  use physics_mod
  use filter_mod
  use test_forcing_mod
  use perf_mod

  implicit none

  private

  public gmcore_init
  public gmcore_init_stage1
  public gmcore_init_stage2
  public gmcore_init_stage3
  public gmcore_run
  public gmcore_final

  public adv_accum_wind
  public block_type
  public dstate_type
  public dtend_type
  public proc

  procedure(space_operators_interface), pointer :: operators

contains

  subroutine gmcore_init(namelist_path, comm)

    character(*), intent(in) :: namelist_path
    integer, intent(in), optional :: comm

    call gmcore_init_stage1(namelist_path, comm)
    call gmcore_init_stage2(namelist_path)

  end subroutine gmcore_init

  subroutine gmcore_init_stage1(namelist_path, comm)

    character(*), intent(in) :: namelist_path
    integer, intent(in), optional :: comm

    call log_init()
    call gas_mixture_init(planet)
    call const_init(planet)
    call time_scheme_init()
    call time_init(dt_dyn)
    call global_mesh%init_global(nlon, nlat, nlev, lon_hw=lon_hw, lat_hw=lat_hw)
    call process_init(comm)
    call process_create_blocks()
    call damp_init()
    call history_init_stage1()

    associate (mesh => blocks(1)%mesh)
    min_lon = mesh%full_lon_deg(mesh%full_ims)
    max_lon = mesh%full_lon_deg(mesh%full_ime)
    min_lat = mesh%full_lat_deg(max(mesh%full_jms, 1))
    max_lat = mesh%full_lat_deg(min(mesh%full_jme, global_mesh%full_nlat))
    end associate

  end subroutine gmcore_init_stage1

  subroutine gmcore_init_stage2(namelist_path)

    character(*), intent(in) :: namelist_path

    character(10) time_value, time_units
    integer iblk
    real(r8) seconds

    if (proc%is_root()) then
      print *, ''
      print *, ' ▄▄▄▄▄▄▄ ▄▄   ▄▄ ▄▄▄▄▄▄▄ ▄▄▄▄▄▄▄ ▄▄▄▄▄▄   ▄▄▄▄▄▄▄ '
      print *, '█       █  █▄█  █       █       █   ▄  █ █       █'
      print *, '█   ▄▄▄▄█       █       █   ▄   █  █ █ █ █    ▄▄▄█'
      print *, '█  █  ▄▄█       █     ▄▄█  █ █  █   █▄▄█▄█   █▄▄▄ '
      print *, '█  █ █  █       █    █  █  █▄█  █    ▄▄  █    ▄▄▄█'
      print *, '█  █▄▄█ █ ██▄██ █    █▄▄█       █   █  █ █   █▄▄▄ '
      print *, '█▄▄▄▄▄▄▄█▄█   █▄█▄▄▄▄▄▄▄█▄▄▄▄▄▄▄█▄▄▄█  █▄█▄▄▄▄▄▄▄█'
      print *, ''
    end if

    call vert_coord_init(namelist_path)
    call restart_init()
    call tracer_init()
    call pgf_init()
    call interp_init()
    call operators_init()
    call physics_init_stage1(namelist_path)
    if (baroclinic .and. physics_suite /= 'cam') call tracer_add_moist()
    call tracer_allocate()
    call adv_init()
    call history_init_stage2()

    operators => space_operators

    time_value = split_string(print_interval, ' ', 1)
    time_units = split_string(print_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days', 'sol')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid print interval ' // trim(print_interval) // '!')
    end select

    call time_add_alert('print', seconds=seconds)

    if (proc%is_root()) call print_namelist()

    do iblk = 1, size(blocks)
      blocks(iblk)%mesh%full_lev  = global_mesh%full_lev
      blocks(iblk)%mesh%half_lev  = global_mesh%half_lev
      blocks(iblk)%mesh%full_dlev = global_mesh%full_dlev
      blocks(iblk)%mesh%half_dlev = global_mesh%half_dlev
      call blocks(iblk)%static%init_stage2(blocks(iblk)%mesh)
    end do

  end subroutine gmcore_init_stage2

  subroutine gmcore_init_stage3()

    call physics_init_stage2()

  end subroutine gmcore_init_stage3

  subroutine gmcore_run()

    integer i, j, m, iblk, itime

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)     , &
                 mesh  => blocks(iblk)%mesh, &
                 dstate => blocks(iblk)%dstate(old))
      if (baroclinic) then
        ! Ensure bottom gz_lev is the same as gzs.
        do itime = lbound(block%dstate, 1), ubound(block%dstate, 1)
          do j = mesh%full_jms, mesh%full_jme
            do i = mesh%full_ims, mesh%full_ime
              block%dstate(itime)%gz_lev%d(i,j,global_mesh%half_kde) = block%static%gzs%d(i,j)
            end do
          end do
        end do
      end if
      call dstate%c2a()
      call calc_div(block, dstate)
      call tracer_calc_qm(block)
      end associate
    end do

    call pgf_init_after_ic()
    call operators_prepare(blocks, old, dt_dyn)
    call adv_prepare(old)
    call diagnose(blocks, old)
    if (proc%is_root()) call log_print_diag(curr_time%isoformat())
    call output(old)

    model_main_loop: do while (.not. time_is_finished())
      ! ------------------------------------------------------------------------
      !                              Dynamical Core
      do iblk = 1, size(blocks)
        call time_integrator(operators, blocks(iblk), old, new, dt_dyn)
        call damp_run(blocks(iblk), blocks(iblk)%dstate(new), dt_dyn)
        if (pdc_type == 1) call physics_update_dynamics(blocks(iblk), new, dt_dyn)
        call blocks(iblk)%dstate(new)%c2a()
      end do
      ! ------------------------------------------------------------------------
      ! Advance to n+1 time level.
      ! NOTE: Time indices are swapped, e.g. new <=> old.
      call time_advance(dt_dyn)
      ! ------------------------------------------------------------------------
      !                            Tracer Advection
      call adv_accum_wind(new)
      call adv_run(new)
      ! ------------------------------------------------------------------------
      !                                Physics
      call test_forcing_run(dt_dyn, new)
      if (baroclinic) then
        do iblk = 1, size(blocks)
          call physics_run(blocks(iblk), new, dt_phys)
          if (pdc_type == 3) call physics_update(blocks(iblk), new, dt_phys)
        end do
      end if
      ! ------------------------------------------------------------------------
      call diagnose(blocks, new)
      if (proc%is_root() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call blocks(1)%accum(new)
      call output(new)
    end do model_main_loop

  end subroutine gmcore_run

  subroutine gmcore_final()

    call log_final()
    call time_final()
    call interp_final()
    call gas_mixture_final()
    call vert_coord_final()
    call tracer_final()
    call pgf_final()
    call adv_final()
    call damp_final()
    call history_final()
    call process_final()
    call perf_final()

  end subroutine gmcore_final

  subroutine output(itime)

    integer, intent(in) :: itime

    logical, save :: first_call = .true.
    real(8), save :: time1, time2

    if (first_call .or. time_is_alerted('history_write')) then
      if (first_call) time1 = MPI_WTIME()
      call process_barrier()
      time2 = MPI_WTIME()
      if (.not. first_call) then
        if (proc%is_root()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      if (output_h0) call history_write_h0(itime)
      if (output_h1) call history_write_h1(itime)
    end if
    if (time_is_alerted('restart_write')) then
      call restart_write(itime)
    end if
    first_call = .false.

  end subroutine output

  subroutine diagnose(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    integer i, j, k, iblk
    real(r8) tm, te, tpe, tpt, max_w
    real(r8) te_ke, te_ie, te_pe

    tm    = 0
    te    = 0
    tpe   = 0
    tpt   = 0
    te_ke = 0
    te_ie = 0
    te_pe = 0
    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                    , &
                 mesh    => blocks(iblk)%mesh               , &
                 gzs     => blocks(iblk)%static%gzs         , &
                 mgs     => blocks(iblk)%dstate(itime)%mgs  , &
                 dmg     => blocks(iblk)%dstate(itime)%dmg  , &
                 dmg_lon => blocks(iblk)%aux%dmg_lon        , &
                 dmg_lat => blocks(iblk)%aux%dmg_lat        , &
                 u_lon  => blocks(iblk)%dstate(itime)%u_lon , &
                 v_lat  => blocks(iblk)%dstate(itime)%v_lat , &
                 tv     => blocks(iblk)%dstate(itime)%tv    , &
                 pt     => blocks(iblk)%dstate(itime)%pt    , &
                 pv_lon => blocks(iblk)%aux%pv_lon          , &
                 pv_lat => blocks(iblk)%aux%pv_lat          )
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            tm = tm + dmg%d(i,j,k) * mesh%area_cell(j)
          end do
        end do
      end do

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            te_ke = te_ke + dmg_lon%d(i,j,k) * 0.5_r8 * u_lon%d(i,j,k)**2 * mesh%area_lon(j) * 2
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            te_ke = te_ke + dmg_lat%d(i,j,k) * 0.5_r8 * v_lat%d(i,j,k)**2 * mesh%area_lat(j) * 2
          end do
        end do
      end do
      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              te_ie = te_ie + dmg%d(i,j,k) * cpd * tv%d(i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            te_pe = te_pe + gzs%d(i,j) * mgs%d(i,j) * mesh%area_cell(j)
          end do
        end do
      else
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            te_pe = te_pe + (dmg%d(i,j,1)**2 * g * 0.5_r8 + dmg%d(i,j,1) * gzs%d(i,j)) * mesh%area_cell(j)
          end do
        end do
      end if

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            tpe = tpe + dmg_lon%d(i,j,k) * pv_lon%d(i,j,k)**2 * 0.5_r8 * mesh%area_lon(j) * 2
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            tpe = tpe + dmg_lat%d(i,j,k) * pv_lat%d(i,j,k)**2 * 0.5_r8 * mesh%area_lat(j) * 2
          end do
        end do
      end do

      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tpt = tpt + dmg%d(i,j,k) * pt%d(i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
      end if
      end associate
    end do
    call global_sum(proc%comm, tm)
    call global_sum(proc%comm, te_ke)
    call global_sum(proc%comm, te_ie)
    call global_sum(proc%comm, te_pe)
    call global_sum(proc%comm, tpe)
    if (baroclinic) call global_sum(proc%comm, tpt)
    te = te_ke + te_ie + te_pe

    do iblk = 1, size(blocks)
      blocks(iblk)%dstate(itime)%tm  = tm
      blocks(iblk)%dstate(itime)%te  = te
      blocks(iblk)%dstate(itime)%tpe = tpe
      blocks(iblk)%dstate(itime)%te_ke = te_ke
      blocks(iblk)%dstate(itime)%te_ie = te_ie
      blocks(iblk)%dstate(itime)%te_pe = te_pe
      ! call calc_omg(blocks(iblk), blocks(iblk)%dstate(itime))
    end do

    call log_add_diag('tm' , tm )
    if (baroclinic) call log_add_diag('tpt', tpt)
    call log_add_diag('te' , te )
    call log_add_diag('tpe', tpe)

  end subroutine diagnose

  subroutine space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, pass, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type), intent(inout) :: tend1
    type(dtend_type), intent(in) :: tend2
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: substep

    integer i, j, k

    call tend1%reset_flags()

    associate (mesh => block%mesh)
    select case (pass)
    case (all_pass)
      call operators_prepare(block, star_state, dt, pass, substep)
      if (hydrostatic) then
        call calc_grad_mf          (block, star_state, tend1, dt)
        call calc_dmgsdt           (block, star_state, tend1, dt)
        call calc_we_lev           (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_grad_ptf         (block, star_state, tend1, dt)
        call calc_coriolis         (block, star_state, tend1, dt)
        call calc_grad_ke          (block, star_state, tend1, dt)
        call pgf_run               (block, star_state, tend1)

        tend1%update_u   = .true.
        tend1%update_v   = .true.
        tend1%update_mgs = .true.
        tend1%update_pt  = .true.
      else
        call calc_grad_mf        (block, star_state, tend1, dt)
        call calc_coriolis       (block, star_state, tend1, dt)
        call calc_grad_ke        (block, star_state, tend1, dt)
        call pgf_run             (block, star_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dgz%d(i,j,k) = -block%aux%dmf%d(i,j,k) * g
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
        tend1%update_gz = .true.
      end if
    case (forward_pass)
      call operators_prepare(block, star_state, dt, pass, substep)
      if (hydrostatic) then
        call calc_grad_mf          (block, star_state, tend1, dt)
        call calc_dmgsdt           (block, star_state, tend1, dt)
        call calc_we_lev           (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_grad_ptf         (block, star_state, tend1, dt)
        call calc_coriolis         (block, star_state, tend1, dt)
        call calc_grad_ke          (block, star_state, tend1, dt)

        tend1%update_mgs = .true.
        tend1%update_pt  = .true.
      else
        call calc_grad_mf         (block, star_state, tend1, dt)
        call calc_coriolis        (block, star_state, tend1, dt)
        call calc_grad_ke         (block, star_state, tend1, dt)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dgz%d(i,j,k) = -block%aux%dmf%d(i,j,k) * g
            end do
          end do
        end do

        tend1%update_gz = .true.
      end if
    case (backward_pass)
      call operators_prepare(block, new_state, dt, pass, substep)
      if (hydrostatic) then
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du%d(i,j,k) = tend1%du%d(i,j,k) + tend2%du%d(i,j,k)
#ifdef OUTPUT_H1_DTEND
              tend1%dudt_coriolis%d(i,j,k) = tend2%dudt_coriolis%d(i,j,k)
              tend1%dudt_wedudeta%d(i,j,k) = tend2%dudt_wedudeta%d(i,j,k)
              tend1%dudt_dkedx   %d(i,j,k) = tend2%dudt_dkedx   %d(i,j,k)
#endif
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv%d(i,j,k) = tend1%dv%d(i,j,k) + tend2%dv%d(i,j,k)
#ifdef OUTPUT_H1_DTEND
              tend1%dvdt_coriolis%d(i,j,k) = tend2%dvdt_coriolis%d(i,j,k)
              tend1%dvdt_wedvdeta%d(i,j,k) = tend2%dvdt_wedvdeta%d(i,j,k)
              tend1%dvdt_dkedy   %d(i,j,k) = tend2%dvdt_dkedy   %d(i,j,k)
#endif
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du%d(i,j,k) = tend1%du%d(i,j,k) + tend2%du%d(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv%d(i,j,k) = tend1%dv%d(i,j,k) + tend2%dv%d(i,j,k)
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
      end if
    end select
    end associate

  end subroutine space_operators

end module gmcore_mod
