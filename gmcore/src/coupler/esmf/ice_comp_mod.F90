module ice_comp_mod

  use esmf
  use comp_wrapper_mod

  implicit none

  private

  public ice_comp_set_services
  public comp_type

contains

  subroutine ice_comp_set_services(comp, rc)

    type(ESMF_GridComp), intent(inout) :: comp
    integer, intent(out) :: rc

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=ice_comp_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=ice_comp_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=ice_comp_final, rc=rc)

    rc = ESMF_SUCCESS

  end subroutine ice_comp_set_services

  subroutine ice_comp_init(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ice_comp_init

  subroutine ice_comp_run(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ice_comp_run

  subroutine ice_comp_final(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ice_comp_final

end module ice_comp_mod
