module stepon

  use shr_kind_mod  , only: r8 => shr_kind_r8
  use dyn_comp      , only: dyn_import_t, dyn_export_t
  use physics_types , only: physics_state, physics_tend
  use physics_buffer, only: physics_buffer_desc
  use ppgrid        , only: begchunk, endchunk
  use camsrfexch    , only: cam_out_t, cam_in_t
  use time_manager  , only: is_first_step, get_step_size

  implicit none

  private

  public stepon_init
  public stepon_run1
  public stepon_run2

contains

  subroutine stepon_init(dyn_in, dyn_out)

    type(dyn_import_t), intent(inout) :: dyn_in
    type(dyn_export_t), intent(inout) :: dyn_out

  end subroutine stepon_init

  subroutine stepon_run1(dtime, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out)

    real(r8), intent(out) :: dtime
    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(dyn_import_t), intent(inout) :: dyn_in
    type(dyn_export_t), intent(inout) :: dyn_out

    dtime = get_step_size()

    ! Move data into phys_state structure.

  end subroutine stepon_run1

  subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
    type(dyn_import_t), intent(inout) :: dyn_in
    type(dyn_export_t), intent(inout) :: dyn_out

    ! Move data from physics to dynamics.

  end subroutine stepon_run2

end module stepon