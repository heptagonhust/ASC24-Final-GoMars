module gmcore_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use process_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use restart_mod
  use block_mod
  use vert_coord_mod
  use time_schemes_mod
  use operators_mod
  use moist_mod
  use interp_mod
  use debug_mod
  use adv_mod
  use pgf_mod
  use damp_mod
  use diag_state_mod
  use physics_mod
  use filter_mod
  use test_forcing_mod

  implicit none

  private

  public gmcore_init
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

    character(10) time_value, time_units
    integer iblk
    real(r8) seconds

    call log_init()
    call global_mesh%init_global(nlon, nlat, nlev, lon_hw=lon_hw, lat_hw=2)
    call process_init(comm)
    call vert_coord_init(nlev, namelist_path)
    call process_create_blocks()
    call time_init(dt_dyn)
    call diag_state_init(blocks)
    call restart_init()
    call time_scheme_init()
    call adv_init()
    call pgf_init()
    call interp_init()
    call operators_init()
    call physics_init()
    call damp_init(blocks)
    if (baroclinic) call moist_init()
    call adv_allocate_tracers(blocks)
    call history_init()

    operators => space_operators

    time_value = split_string(print_interval, ' ', 1)
    time_units = split_string(print_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
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

    if (is_root_proc()) call print_namelist()

    do iblk = 1, size(blocks)
      if (baroclinic) call moist_link_state(blocks(iblk))
    end do

  end subroutine gmcore_init

  subroutine gmcore_run()

    integer m, iblk, itime

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)     , &
                 mesh  => blocks(iblk)%mesh, &
                 dstate => blocks(iblk)%dstate(old))
      if (baroclinic) then 
        call prepare_static(block)
        ! Ensure bottom gz_lev is the same as gzs.
        do itime = lbound(block%dstate, 1), ubound(block%dstate, 1)
          block%dstate(itime)%gz_lev(:,:,global_mesh%half_kde) = block%static%gzs
        end do
      end if
      call blocks(iblk)%dstate(old)%c2a()
      if (baroclinic) call moist_link_state(block)
      end associate
    end do

    call operators_prepare(blocks, old, dt_dyn)
    call adv_prepare(old)
    if (nonhydrostatic) call nh_prepare(blocks)
    call diagnose(blocks, old)
    call output(old)
    if (is_root_proc()) call log_print_diag(curr_time%isoformat())

    model_main_loop: do while (.not. time_is_finished())
      ! ------------------------------------------------------------------------
      !                              Dynamical Core
      do iblk = 1, size(blocks)
        call time_integrator(operators, blocks(iblk), old, new, dt_dyn)
        call damp_run(blocks(iblk), blocks(iblk)%dstate(new), blocks(iblk)%dtend(new), dt_dyn)
        call blocks(iblk)%dstate(new)%c2a()
      end do

      ! Advance to n+1 time level.
      ! NOTE: Time indices are swapped, e.g. new <=> old.
      call time_advance(dt_dyn)
      ! ------------------------------------------------------------------------
      !                            Tracer Advection
      do iblk = 1, size(blocks)
        call adv_run(blocks(iblk), old)
      end do
      ! ------------------------------------------------------------------------
      !                                Physics
      if (baroclinic) then
        do iblk = 1, size(blocks)
          call test_forcing_run(blocks(iblk), dt_dyn, blocks(iblk)%static, blocks(iblk)%dstate(old))
          call moist_link_state(blocks(iblk))
          call physics_run_after_dynamics(blocks(iblk), old, dt_phys)
        end do
      end if
      ! ------------------------------------------------------------------------
      call diagnose(blocks, old)
      if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call output(old)
    end do model_main_loop

    ! Write a restart file at last.
    ! call restart_write(old)

  end subroutine gmcore_run

  subroutine gmcore_final()

    call log_final()
    call time_final()
    call interp_final()
    call adv_final()
    call damp_final()
    call diag_state_final()
    call history_final()
    call process_final()

  end subroutine gmcore_final

  subroutine prepare_static(block)

    class(block_type), intent(inout) :: block

    integer i, j

    associate (mesh    => block%mesh          , &
               gzs     => block%static%gzs    , & ! in
               dzsdlon => block%static%dzsdlon, & ! out
               dzsdlat => block%static%dzsdlat)   ! out
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          dzsdlon(i,j) = (gzs(i+1,j) - gzs(i,j)) / g / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          dzsdlat(i,j) = (gzs(i,j+1) - gzs(i,j)) / g / mesh%de_lat(j)
        end do
      end do
      call fill_halo(block%halo, dzsdlon, full_lon=.false., full_lat=.true.)
      call fill_halo(block%halo, dzsdlat, full_lon=.true., full_lat=.false.)
    end associate

  end subroutine prepare_static

  subroutine output(itime)

    integer, intent(in) :: itime

    real(r8), save :: time1 = 0, time2
    integer i, j, k, iblk

    if (time_step == 0 .or. time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      if (output_h0) call history_write_h0(blocks, itime)
      if (output_h1) call history_write_h1(blocks, itime)
    end if
    if (time_is_alerted('restart_write')) then
      call restart_write(itime)
    end if

  end subroutine output

  subroutine diagnose(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    integer i, j, k, iblk
    real(r8) tm, te, tav, tpe, tpt, max_w
    real(r8) te_ke, te_ie, te_pe

    tm    = 0
    te    = 0
    tav   = 0
    tpe   = 0
    tpt   = 0
    te_ke = 0
    te_ie = 0
    te_pe = 0
    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh, &
                 dstate  => blocks(iblk)%dstate(itime), &
                 static => blocks(iblk)%static)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            tm = tm + dstate%m(i,j,k) * mesh%area_cell(j)
          end do
        end do
      end do

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            te_ke = te_ke + dstate%mfx_lon(i,j,k) * 0.5_r8 * dstate%u_lon(i,j,k) * mesh%area_lon(j) * 2
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            te_ke = te_ke + dstate%mfy_lat(i,j,k) * 0.5_r8 * dstate%v_lat(i,j,k) * mesh%area_lat(j) * 2
          end do
        end do
      end do
      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              te_ie = te_ie + dstate%m(i,j,k) * cpd * dstate%t(i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            te_pe = te_pe + static%gzs(i,j) * dstate%phs(i,j) * mesh%area_cell(j)
          end do
        end do
      else
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            te_pe = te_pe + (dstate%m(i,j,1)**2 * g * 0.5_r8 + dstate%m(i,j,1) * static%gzs(i,j)) * mesh%area_cell(j)
          end do
        end do
      end if

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids, mesh%half_ide
            tav = tav + dstate%pv(i,j,k) * mesh%area_vtx(j)
          end do
        end do
      end do

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            tpe = tpe + dstate%m_lon(i,j,k) * dstate%pv_lon(i,j,k)**2 * 0.5_r8 * mesh%area_lon(j)
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            tpe = tpe + dstate%m_lat(i,j,k) * dstate%pv_lat(i,j,k)**2 * 0.5_r8 * mesh%area_lat(j)
          end do
        end do
      end do

      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tpt = tpt + dstate%m(i,j,k) * dstate%pt(i,j,k) * mesh%area_cell(j)
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
    call global_sum(proc%comm, tav)
    call global_sum(proc%comm, tpe)
    if (baroclinic) call global_sum(proc%comm, tpt)
    te = te_ke + te_ie + te_pe

    do iblk = 1, size(blocks)
      blocks(iblk)%dstate(itime)%tm  = tm
      blocks(iblk)%dstate(itime)%te  = te
      blocks(iblk)%dstate(itime)%tav = tav
      blocks(iblk)%dstate(itime)%tpe = tpe
      blocks(iblk)%dstate(itime)%te_ke = te_ke
      blocks(iblk)%dstate(itime)%te_ie = te_ie
      blocks(iblk)%dstate(itime)%te_pe = te_pe
      if (diag_state(iblk)%is_init()) call diag_state(iblk)%run(blocks(iblk)%dstate(itime))
    end do

    call log_add_diag('tm' , tm )
    if (baroclinic) call log_add_diag('tpt', tpt)
    call log_add_diag('te' , te )
    call log_add_diag('tpe', tpe)

    if (nonhydrostatic) then
      max_w = 0
      do iblk = 1, size(blocks)
        max_w = max(max_w, maxval(abs(blocks(iblk)%dstate(itime)%w)))
      end do
      call global_max(proc%comm, max_w)
      call log_add_diag('w', max_w)
    end if

  end subroutine diagnose

  subroutine space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, pass)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    type(dtend_type), intent(inout) :: tend1
    type(dtend_type), intent(in) :: tend2
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    integer i, j, k

    call tend1%reset_flags()

    associate (mesh => block%mesh)
    select case (pass)
    case (all_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_grad_mf          (block, star_state, tend1, dt)
        call calc_dphsdt           (block, star_state, tend1, dt)
        call calc_we_lev           (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_grad_ptf         (block, star_state, tend1, dt)
        call calc_coriolis         (block, star_state, tend1, dt)
        call calc_grad_ke          (block, star_state, tend1, dt)
        call pgf_run               (block, star_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then
        call calc_grad_mf          (block, star_state, tend1, dt)
        call calc_dphsdt           (block, star_state, tend1, dt)
        call calc_we_lev           (block, star_state, tend1, dt)
        call calc_grad_ptf         (block, star_state, tend1, dt)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
        call update_state(block, tend1, old_state, new_state, dt)

        call nh_solve(block, tend1, old_state, star_state, new_state, dt)

        call calc_coriolis         (block, star_state, tend1, dt)
        call calc_grad_ke          (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call pgf_run               (block,  new_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else
        call calc_grad_mf        (block, star_state, tend1, dt)
        call calc_coriolis       (block, star_state, tend1, dt)
        call calc_grad_ke        (block, star_state, tend1, dt)
        call pgf_run             (block, star_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dgz(i,j,k) = - (tend1%dmfdlon(i,j,k) + tend1%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
        tend1%update_gz = .true.
      end if
    case (forward_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_grad_mf          (block, star_state, tend1, dt)
        call calc_dphsdt           (block, star_state, tend1, dt)
        call calc_we_lev           (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_grad_ptf         (block, star_state, tend1, dt)
        call calc_coriolis         (block, star_state, tend1, dt)
        call calc_grad_ke          (block, star_state, tend1, dt)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then

      else
        call calc_grad_mf         (block, star_state, tend1, dt)
        call calc_coriolis        (block, star_state, tend1, dt)
        call calc_grad_ke         (block, star_state, tend1, dt)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dgz(i,j,k) = - (tend1%dmfdlon(i,j,k) + tend1%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend1%update_gz = .true.
      end if
    case (backward_pass)
      call operators_prepare(block, new_state, dt, pass)
      if (hydrostatic) then
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) = tend2%du(i,j,k) - tend1%pgf_lon(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = tend2%dv(i,j,k) - tend1%pgf_lat(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else if (nonhydrostatic) then

      else
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids, mesh%half_ide
              tend1%du(i,j,k) = tend2%du(i,j,k) - tend1%pgf_lon(i,j,k)
            end do
          end do

          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              tend1%dv(i,j,k) = tend2%dv(i,j,k) - tend1%pgf_lat(i,j,k)
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
