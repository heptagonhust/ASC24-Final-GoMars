module comp_wrapper_mod

  use esmf
  use string

  implicit none

  private

  public comp_info_type
  public comp_info_wrapper_type
  public comp_type
  public set_services_interface
  public add_field

  type comp_info_type
  sequence
    character(256) :: namelist_path = 'N/A'
    type(ESMF_Grid) grid
  end type comp_info_type

  type comp_info_wrapper_type
  sequence
    type(comp_info_type), pointer :: p
  end type comp_info_wrapper_type

  type comp_type
    character(30) :: name = 'N/A'
    type(ESMF_GridComp) comp
    type(ESMF_State) import_state
    type(ESMF_State) export_state
    type(ESMF_Clock) clock
    type(comp_info_type) info
  contains
    procedure :: init => comp_init
    procedure :: run => comp_run
    procedure :: clear => comp_clear
    final :: comp_final
  end type comp_type

  interface
    subroutine set_services_interface(grid_comp, rc)
      import ESMF_GridComp
      type(ESMF_GridComp) grid_comp
      integer, intent(out) :: rc
    end subroutine set_services_interface
  end interface

contains

  subroutine comp_init(this, comp_name, set_services, clock, namelist_path, phase)

    class(comp_type), intent(inout), target :: this
    character(*), intent(in) :: comp_name
    procedure(set_services_interface) set_services
    type(ESMF_Clock), intent(inout) :: clock
    character(*), intent(in) :: namelist_path
    integer, intent(in), optional :: phase

    type(comp_info_wrapper_type) info_wrapper
    integer rc

    this%name = comp_name
    this%clock = clock
    this%info%namelist_path = namelist_path

    this%comp = ESMF_GridCompCreate(contextflag=ESMF_CONTEXT_PARENT_VM, name=comp_name, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Component ' // trim(comp_name) // ' is created.', ESMF_LOGMSG_INFO)

    call ESMF_GridCompSetServices(this%comp, userRoutine=set_services, rc=rc)
    if (ESMF_LogFoundError(rc, msg='Registration failed', rcToReturn=rc)) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    info_wrapper%p => this%info
    call ESMF_GridCompSetInternalState(this%comp, info_wrapper, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    this%import_state = ESMF_StateCreate(name=comp_name, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Import state of ' // trim(comp_name) // ' is created.', ESMF_LOGMSG_INFO)

    this%export_state = ESMF_StateCreate(name=comp_name, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Export state of ' // trim(comp_name) // ' is created.', ESMF_LOGMSG_INFO)

    call ESMF_LogWrite('Start to initialize component ' // trim(comp_name) // '.', ESMF_LOGMSG_INFO)
    call ESMF_GridCompInitialize(this%comp, importState=this%import_state, exportState=this%export_state, clock=clock, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

  end subroutine comp_init

  subroutine comp_run(this, phase)

    class(comp_type), intent(inout) :: this
    integer, intent(in), optional :: phase

    integer p, rc

    p = merge(phase, 1, present(phase))

    call ESMF_LogWrite('Run component ' // trim(this%name) // ' at phase ' // to_str(p) // '.', ESMF_LOGMSG_INFO)
    call ESMF_GridCompRun(this%comp, importState=this%import_state, exportState=this%export_state, clock=this%clock, phase=p, rc=rc)

  end subroutine comp_run

  subroutine comp_clear(this)

    class(comp_type), intent(inout) :: this

    integer rc

    call ESMF_GridCompDestroy(this%comp, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Component ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

    call ESMF_GridDestroy(this%info%grid, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Grid of ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(this%import_state, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Import state of ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(this%export_state, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Export state of ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

  end subroutine comp_clear

  subroutine comp_final(this)

    type(comp_type), intent(inout) :: this

    call this%clear()

  end subroutine comp_final

  subroutine add_field(info, state, field_name)

    type(comp_info_wrapper_type), intent(in) :: info
    type(ESMF_State), intent(inout) :: state
    character(*), intent(in) :: field_name

    type(ESMF_Field) field
    integer rc

    field = ESMF_FieldCreate( &
      grid=info%p%grid, &
      typekind=ESMF_TYPEKIND_R8, &
      indexflag=ESMF_INDEX_GLOBAL, &
      staggerloc=ESMF_STAGGERLOC_CENTER, &
      name=field_name, &
      rc=rc)
    call ESMF_StateAdd(state, [field], rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Add field ' // trim(field_name) // '.', ESMF_LOGMSG_INFO)

  end subroutine add_field

end module comp_wrapper_mod
