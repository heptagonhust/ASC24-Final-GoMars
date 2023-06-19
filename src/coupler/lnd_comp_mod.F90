module lnd_comp_mod

  use esmf
  use comp_wrapper_mod

  implicit none

  private

  public lnd_comp_SetServices
  public comp_type

contains

  subroutine lnd_comp_SetServices(comp, rc)

    type(ESMF_GridComp) comp
    integer, intent(out) :: rc

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=lnd_comp_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=lnd_comp_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=lnd_comp_final, rc=rc)

    rc = ESMF_SUCCESS

  end subroutine lnd_comp_SetServices

  subroutine lnd_comp_init(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine lnd_comp_init

  subroutine lnd_comp_run(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine lnd_comp_run

  subroutine lnd_comp_final(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine lnd_comp_final

end module lnd_comp_mod
