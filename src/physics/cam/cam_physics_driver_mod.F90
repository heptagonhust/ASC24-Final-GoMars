module cam_physics_driver_mod

  use shr_kind_mod    , only: r8 => shr_kind_r8
  use spmd_utils      , only: spmdinit
  use check_energy    , only: check_energy_timestep_init
  use dyn_grid        , only: dyn_grid_init, dyn_grid_final
  use dyn_comp        , only: dyn_init, dyn_final, dyn_import_t, dyn_export_t
  use phys_grid       , only: phys_grid_init, local_dp_map, get_lon_all_p, get_lat_all_p
  use physics_types   , only: physics_state, physics_tend, set_state_pdry
  use physics_buffer  , only: physics_buffer_desc, pbuf_get_chunk
  use physpkg         , only: phys_register, phys_init, phys_run1, phys_run2, phys_final
  use physconst       , only: cappa, cpair, gravit
  use ppgrid          , only: begchunk, endchunk, pcols, pver, pverp
  use qneg_module     , only: qneg3
  use constituents    , only: pcnst, cnst_name, cnst_longname, cnst_get_type_byind, qmin
  use chem_surfvals   , only: chem_surfvals_init
  use air_composition , only: air_composition_init
  use camsrfexch      , only: cam_out_t, cam_in_t, hub2atm_alloc, atm2hub_alloc
  use runtime_opts    , only: read_namelist
  use shr_cal_mod     , only: shr_cal_gregorian
  use shr_pio_mod     , only: shr_pio_init1, shr_pio_init2
  use cam_pio_utils   , only: init_pio_subsystem
  use cam_instance    , only: cam_instance_init
  use cam_initfiles   , only: cam_initfiles_open
  use cam_control_mod , only: cam_ctrl_init
  use time_manager    , only: timemgr_init

  use flogger
  use string          , only: to_int
  use namelist_mod    , only: dt_phys, dt_adv, restart, cam_namelist_path, case_name, case_desc
  use process_mod     , only: proc, process_barrier
  use time_mod        , only: start_time, end_time, curr_time
  use tracer_mod      , only: tracer_add, tracer_get_array_qm, tracers
  use block_mod       , only: block_type

  implicit none

  private

  public cam_physics_driver_init
  public cam_physics_driver_final
  public cam_physics_driver_run1
  public cam_physics_driver_run2
  public cam_physics_d2p
  public cam_physics_p2d

  type(physics_state), pointer :: phys_state(:) => null()
  type(physics_tend), pointer :: phys_tend(:) => null()
  type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

  type(dyn_import_t) dyn_in
  type(dyn_export_t) dyn_out
  type(cam_in_t), pointer :: cam_in(:) => null()
  type(cam_out_t), pointer :: cam_out(:) => null()

  integer, parameter :: atm_id = 1
  integer, parameter :: ncomps = 1
  integer, parameter :: comp_id(ncomps) = [atm_id]
  character(3), parameter :: comp_name(ncomps) = ['atm']
  logical comp_iamin(ncomps)
  integer comp_comm(ncomps)
  integer comp_comm_iam(ncomps)
  logical, allocatable :: wetq(:)

contains

  subroutine cam_physics_driver_init(namelist_path)

    character(*), intent(in) :: namelist_path

    integer i

    comp_comm(1) = proc%comm
    comp_comm_iam(1) = proc%id
    comp_iamin(1) = .true.

    call spmdinit(proc%comm)
    call cam_instance_init(atm_id, 'atm', 1, '_1')
    call shr_pio_init1(ncomps, merge(namelist_path, cam_namelist_path, cam_namelist_path == 'N/A'), proc%comm)
    call shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
    call init_pio_subsystem()
    call cam_ctrl_init(              &
      caseid_in=case_name          , &
      ctitle_in=case_desc          , &
      initial_run_in=.not.restart  , &
      restart_run_in=restart       , &
      branch_run_in=.false.        , &
      post_assim_in=.false.        , &
      aqua_planet_in=.false.       , &
      brnch_retain_casename_in=.true.)
    call timemgr_init( &
      dtime_in=int(dt_phys), &
      calendar_in=shr_cal_gregorian, &
      start_ymd=to_int(start_time%format('%y%m%d')), &
      start_tod=start_time%hour*3600+start_time%minute*60+start_time%second, &
      ref_ymd=20000101, &
      ref_tod=0, &
      stop_ymd=to_int(end_time%format('%y%m%d')), &
      stop_tod=end_time%hour*3600+end_time%minute*60+end_time%second, &
      curr_ymd=to_int(curr_time%format('%y%m%d')), &
      curr_tod=curr_time%hour*3600+curr_time%minute*60+curr_time%second, &
      perpetual_run=.false., &
      perpetual_ymd=20000101, &
      initial_run=.not. restart &
    )
    call read_namelist(merge(namelist_path, cam_namelist_path, cam_namelist_path == 'N/A'))
    call cam_initfiles_open()
    call dyn_grid_init()
    call phys_grid_init()
    call phys_register()
    call dyn_init(dyn_in, dyn_out)
    call chem_surfvals_init()
    ! call air_composition_init()

    allocate(wetq(pcnst))
    do i = 1, pcnst
      call tracer_add('cam_cnst', dt_adv, cnst_name(i), cnst_longname(i), 'kg kg-1', type=0)
      wetq(i) = cnst_get_type_byind(i) == 'wet'
    end do

    if (.not. restart) then
      call hub2atm_alloc(cam_in)
      call atm2hub_alloc(cam_out)
    end if

    call phys_init(phys_state, phys_tend, pbuf2d, cam_in, cam_out)

  end subroutine cam_physics_driver_init

  subroutine cam_physics_driver_run1()

    call phys_run1(dt_phys, phys_state, phys_tend, pbuf2d, cam_in, cam_out)

  end subroutine cam_physics_driver_run1

  subroutine cam_physics_driver_run2()

    call phys_run2(phys_state, dt_phys, phys_tend, pbuf2d,  cam_out, cam_in)

  end subroutine cam_physics_driver_run2

  subroutine cam_physics_driver_final()

    call phys_final(phys_state, phys_tend , pbuf2d)
    call dyn_grid_final()
    call dyn_final()

    if (allocated(wetq)) deallocate(wetq)

  end subroutine cam_physics_driver_final

  subroutine cam_physics_d2p(block, itime)

    type(block_type), intent(in) :: block
    integer, intent(in) :: itime

    integer lchnk, ncol, i, k, m
    integer ilon(pcols), jlat(pcols)
    real(r8), pointer, dimension(:,:,:) :: qm
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    call tracer_get_array_qm(block%id, qm)

    associate (mesh => block%mesh, dstate => block%dstate(itime), static => block%static, aux => block%aux)
    if (local_dp_map) then
      do lchnk = begchunk, endchunk
        ncol = phys_state(lchnk)%ncol
        call get_lon_all_p(lchnk, ncol, ilon)
        call get_lat_all_p(lchnk, ncol, jlat)
        do i = 1, ncol
          phys_state(lchnk)%ps(i) = dstate%ps(ilon(i),jlat(i))
          phys_state(lchnk)%phis(i) = static%gzs(ilon(i),jlat(i))
        end do
        do k = 1, mesh%full_nlev
          do i = 1, ncol
            phys_state(lchnk)%u(i,k) = dstate%u(ilon(i),jlat(i),k)
            phys_state(lchnk)%v(i,k) = dstate%v(ilon(i),jlat(i),k)
            phys_state(lchnk)%t(i,k) = dstate%t(ilon(i),jlat(i),k)
            phys_state(lchnk)%exner(i,k) = (dstate%p_lev(ilon(i),jlat(i),mesh%half_nlev) / dstate%p(ilon(i),jlat(i),k))**cappa
            phys_state(lchnk)%omega(i,k) = aux%omg(ilon(i),jlat(i),k)
            phys_state(lchnk)%pmid(i,k) = dstate%p(ilon(i),jlat(i),k)
            phys_state(lchnk)%pdel(i,k) = dstate%p_lev(ilon(i),jlat(i),k+1) - dstate%p_lev(ilon(i),jlat(i),k)
            phys_state(lchnk)%pdeldry(i,k) = dstate%dmg(ilon(i),jlat(i),k)
            phys_state(lchnk)%rpdel(i,k)  = 1.0_r8 / phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
            phys_state(lchnk)%zm(i,k) = dstate%gz(ilon(i),jlat(i),k) / gravit
          end do
        end do
        do k = 1, mesh%half_nlev
          do i = 1, ncol
            phys_state(lchnk)%pint(i,k) = dstate%p_lev(ilon(i),jlat(i),k)
            phys_state(lchnk)%lnpint(i,k) = log(phys_state(lchnk)%pint(i,k))
            phys_state(lchnk)%zi(i,k) = dstate%gz_lev(ilon(i),jlat(i),k) / gravit
          end do
        end do
        do m = 1, pcnst
          if (wetq(m)) then
            do k = 1, mesh%full_nlev
              do i = 1, ncol
                phys_state(lchnk)%q(i,k,m) = tracers(block%id)%q(ilon(i),jlat(i),k,m) / (1 + qm(ilon(i),jlat(i),k))
              end do
            end do
          else
            do k = 1, mesh%full_nlev
              do i = 1, ncol
                phys_state(lchnk)%q(i,k,m) = tracers(block%id)%q(ilon(i),jlat(i),k,m)
              end do
            end do
          end if
        end do
      end do
    else
      call log_error('cam_physics_d2p: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if
    end associate

    do lchnk = begchunk, endchunk
      ncol = phys_state(lchnk)%ncol
      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
        do i = 1, ncol
           phys_state(lchnk)%s(i,k) = cpair * phys_state(lchnk)%t(i,k) &
                                    + gravit * phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
        end do
      end do

      call set_state_pdry(phys_state(lchnk), pdeld_calc=.false.)

      call qneg3('cam_physics_d2p', lchnk, ncol, pcols, pver, 1, pcnst, qmin, phys_state(lchnk)%q)

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)
    end do

  end subroutine cam_physics_d2p

  subroutine cam_physics_p2d(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer lchnk, ncol, i, k, m
    integer ilon(pcols), jlat(pcols)

    associate (ptend => block%ptend)
    ptend%updated_u = .true.
    ptend%updated_v = .true.
    ptend%updated_t = .true.
    if (local_dp_map) then
      do lchnk = begchunk, endchunk
        ncol = phys_state(lchnk)%ncol
        call get_lon_all_p(lchnk, ncol, ilon)
        call get_lat_all_p(lchnk, ncol, jlat)
        do k = 1, pver
          do i = 1, ncol
            ptend%dtdt(i,k) = phys_tend(lchnk)%dtdt(i,k)
            ptend%dudt(i,k) = phys_tend(lchnk)%dudt(i,k)
            ptend%dvdt(i,k) = phys_tend(lchnk)%dvdt(i,k)
          end do
        end do
      end do
    else
      call log_error('cam_physics_p2d: Distributed physics columns are not supported yet!', __FILE__, __LINE__)
    end if
    end associate

  end subroutine cam_physics_p2d

end module cam_physics_driver_mod