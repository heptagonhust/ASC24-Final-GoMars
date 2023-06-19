module ocn_comp_mod

  use esmf
  use comp_wrapper_mod

  implicit none

  private

  public ocn_comp_SetServices
  public comp_type

contains

  subroutine ocn_comp_SetServices(comp, rc)

    type(ESMF_GridComp), intent(inout) :: comp
    integer, intent(out) :: rc

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=ocn_comp_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=ocn_comp_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=ocn_comp_final, rc=rc)

    rc = ESMF_SUCCESS

  end subroutine ocn_comp_SetServices

  subroutine ocn_comp_init(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ocn_comp_init

  subroutine ocn_comp_run(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ocn_comp_run

  subroutine ocn_comp_final(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ocn_comp_final

end module ocn_comp_mod
