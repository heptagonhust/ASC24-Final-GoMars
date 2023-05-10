module cam_physics_driver_mod

  use spmd_utils    , only: spmdinit
  use dyn_grid      , only: dyn_grid_init, dyn_grid_final
  use phys_grid     , only: phys_grid_init
  use physics_types , only: physics_state, physics_tend
  use physics_buffer, only: physics_buffer_desc
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

  type(physics_state), pointer :: phys_state(:) => null()
  type(physics_tend), pointer :: phys_tend(:) => null()
  type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

contains

  subroutine cam_physics_driver_init(namelist_path)

    character(*), intent(in) :: namelist_path

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
    call spmdinit(proc%comm)
    call dyn_grid_init()
    call phys_grid_init()

  end subroutine cam_physics_driver_init

  subroutine cam_physics_driver_final()

    call dyn_grid_final()

  end subroutine cam_physics_driver_final

end module cam_physics_driver_mod