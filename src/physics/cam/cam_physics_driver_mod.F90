module cam_physics_driver_mod

  use shr_kind_mod  , only: r8 => shr_kind_r8
  use spmd_utils    , only: spmdinit
  use dyn_grid      , only: dyn_grid_init, dyn_grid_final
  use dyn_comp      , only: dyn_init, dyn_final, dyn_import_t, dyn_export_t
  use phys_grid     , only: phys_grid_init
  use physics_types , only: physics_state, physics_tend
  use physics_buffer, only: physics_buffer_desc
  use physpkg       , only: phys_register, phys_init, phys_run1, phys_run2
  use ppgrid        , only: begchunk, endchunk
  use chem_surfvals , only: chem_surfvals_init
  use camsrfexch    , only: cam_out_t, cam_in_t, hub2atm_alloc, atm2hub_alloc
  use runtime_opts  , only: read_namelist
  use stepon        , only: stepon_init, stepon_run1, stepon_run2
  use shr_cal_mod   , only: shr_cal_gregorian
  use time_manager  , only: timemgr_init

  use string        , only: to_int
  use namelist_mod  , only: dt_phys, restart
  use process_mod   , only: proc
  use time_mod      , only: start_time, end_time, curr_time

  implicit none

  private

  public cam_physics_driver_init
  public cam_physics_driver_final
  public cam_physics_driver_run1

  type(physics_state), pointer :: phys_state(:) => null()
  type(physics_tend), pointer :: phys_tend(:) => null()
  type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

  type(dyn_import_t) dyn_in
  type(dyn_export_t) dyn_out
  type(cam_in_t), pointer :: cam_in(:) => null()
  type(cam_out_t), pointer :: cam_out(:) => null()

  real(r8) dtime_phys ! Time step for physics tendencies. Set by call to stepon_run1, then passed to phys_run*.

contains

  subroutine cam_physics_driver_init(namelist_path)

    character(*), intent(in) :: namelist_path

    call spmdinit(proc%comm)
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
    call read_namelist(namelist_path)
    call dyn_grid_init()
    call phys_grid_init()
    call phys_register()
    call dyn_init(dyn_in, dyn_out)
    call chem_surfvals_init()

    if (.not. restart) then
      call hub2atm_alloc(cam_in)
      call atm2hub_alloc(cam_out)
    end if

    call phys_init(phys_state, phys_tend, pbuf2d, cam_in, cam_out)

    call stepon_init(dyn_in, dyn_out)

  end subroutine cam_physics_driver_init

  subroutine cam_physics_driver_run1()

    call stepon_run1(dtime_phys, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out)

    call phys_run1(dtime_phys, phys_state, phys_tend, pbuf2d, cam_in, cam_out)

  end subroutine cam_physics_driver_run1

  subroutine cam_physics_driver_final()

    call dyn_grid_final()
    call dyn_final()

  end subroutine cam_physics_driver_final

end module cam_physics_driver_mod