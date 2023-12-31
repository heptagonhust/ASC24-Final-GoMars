module lsm_noahmp_driver_mod

  use lsm_noahmp_types_mod
  use NoahmpDriverMainMod

  implicit none

  private

contains

  subroutine lsm_noahmp_run(state)

    type(lsm_noahmp_state_type), intent(inout) :: state

    ! Set julian, yr, month, day (the last two are added by me).

    call NoahmpDriverMain(state)

  end subroutine lsm_noahmp_run

end module lsm_noahmp_driver_mod
