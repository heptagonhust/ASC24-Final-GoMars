module atm_comp_mod

  use mpi
  use esmf

  use namelist_mod
  use comp_wrapper_mod
  use gmcore_mod
  use prepare_mod

  implicit none

  private

  public atm_comp_SetServices
  public comp_type

contains

  subroutine atm_comp_SetServices(comp, rc)

    type(ESMF_GridComp) comp
    integer, intent(out) :: rc

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=atm_comp_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atm_comp_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=atm_comp_final, rc=rc)

    rc = ESMF_SUCCESS

  end subroutine atm_comp_SetServices

  subroutine atm_comp_init(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    type(comp_info_wrapper_type) info

    call ESMF_GridCompGetInternalState(comp, info, rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    call parse_namelist(info%p%namelist_path)
    call gmcore_init_stage1(info%p%namelist_path, comm=MPI_COMM_WORLD)
    call prepare_topo()
    call gmcore_init_stage2(info%p%namelist_path)
    call prepare_bkg()
    call prepare_final()

    ! Add fields to importState.
    ! call ESMF_StateAdd()
    ! Add fields to exportState.

    rc = ESMF_SUCCESS

  end subroutine atm_comp_init

  subroutine atm_comp_run(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine atm_comp_run

  subroutine atm_comp_final(comp, importState, exportState, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) importState
    type(ESMF_State) exportState
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    call gmcore_final()

    rc = ESMF_SUCCESS

  end subroutine atm_comp_final

end module atm_comp_mod
