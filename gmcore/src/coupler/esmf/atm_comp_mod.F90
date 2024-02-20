module atm_comp_mod

  use mpi
  use esmf

  use namelist_mod
  use comp_wrapper_mod
  use gmcore_mod
  use prepare_mod

  implicit none

  private

  public atm_comp_set_services
  public comp_type

contains

  subroutine atm_comp_set_services(comp, rc)

    type(ESMF_GridComp) comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, userRoutine=atm_comp_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atm_comp_run1, phase=1, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, userRoutine=atm_comp_run2, phase=2, rc=rc)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, userRoutine=atm_comp_final, rc=rc)

  end subroutine atm_comp_set_services

  subroutine atm_comp_init(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    type(comp_info_wrapper_type) info

    rc = ESMF_SUCCESS

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

    info%p%grid = ESMF_GridCreate1PeriDim( &
      maxIndex=[global_mesh%full_nlon,global_mesh%full_nlat], &
      regDecomp=[nproc_lon(1),nproc_lat(1)], &
      name='atm_grid', &
      indexflag=ESMF_INDEX_GLOBAL, &
      rc=rc)

    call add_field(info, export_state, 'Sa_topo'   )
    call add_field(info, export_state, 'Sa_z'      )
    call add_field(info, export_state, 'Sa_u'      )
    call add_field(info, export_state, 'Sa_v'      )
    call add_field(info, export_state, 'Sa_tbot'   )
    call add_field(info, export_state, 'Sa_ptem'   )
    call add_field(info, export_state, 'Sa_shum'   )
    call add_field(info, export_state, 'Sa_dens'   )
    call add_field(info, export_state, 'Sa_pbot'   )
    call add_field(info, export_state, 'Sa_pslv'   )
    call add_field(info, export_state, 'Faxa_lwdn' )
    call add_field(info, export_state, 'Faxa_rainc')
    call add_field(info, export_state, 'Faxa_rainl')
    call add_field(info, export_state, 'Faxa_snowc')
    call add_field(info, export_state, 'Faxa_snowl')
    call add_field(info, export_state, 'Faxa_swndr')
    call add_field(info, export_state, 'Faxa_swvdr')
    call add_field(info, export_state, 'Faxa_swndf')
    call add_field(info, export_state, 'Faxa_swvdf')

    call add_field(info, import_state, 'Sx_anidr'  )
    call add_field(info, import_state, 'Sx_anidf'  )
    call add_field(info, import_state, 'Sx_avsdr'  )
    call add_field(info, import_state, 'Sx_avsdf'  )
    call add_field(info, import_state, 'Sl_lfrac'  )
    call add_field(info, import_state, 'Si_ifrac'  )
    call add_field(info, import_state, 'So_ofrac'  )
    call add_field(info, import_state, 'Sx_tref'   )
    call add_field(info, import_state, 'Sx_qref'   )
    call add_field(info, import_state, 'Sx_t'      )
    call add_field(info, import_state, 'So_t'      )
    call add_field(info, import_state, 'Sl_fv'     )
    call add_field(info, import_state, 'Sl_ram1'   )
    call add_field(info, import_state, 'Sl_snowh'  )
    call add_field(info, import_state, 'Si_snowh'  )
    call add_field(info, import_state, 'So_ssq'    )
    call add_field(info, import_state, 'So_re'     )
    call add_field(info, import_state, 'Sx_u10'    )
    call add_field(info, import_state, 'Faxx_taux' )
    call add_field(info, import_state, 'Faxx_tauy' )
    call add_field(info, import_state, 'Faxx_lat'  )
    call add_field(info, import_state, 'Faxx_sen'  )
    call add_field(info, import_state, 'Faxx_lwup' )
    call add_field(info, import_state, 'Faxx_evap' )

  end subroutine atm_comp_init

  subroutine atm_comp_run1(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call gmcore_run()

  end subroutine atm_comp_run1

  subroutine atm_comp_run2(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine atm_comp_run2

  subroutine atm_comp_final(comp, import_state, export_state, clock, rc)

    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call gmcore_final()

  end subroutine atm_comp_final

end module atm_comp_mod
