module comp_wrapper_mod

  use esmf

  implicit none

  private

  public comp_info_type
  public comp_info_wrapper_type
  public comp_type
  public set_services_interface

  type comp_info_type
  sequence
    character(256) :: namelist_path = 'N/A'
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
    type(comp_info_type) info
  contains
    procedure :: init => comp_init
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

  subroutine comp_init(this, comp_name, set_services, clock, namelist_path)

    class(comp_type), intent(inout), target :: this
    character(*), intent(in) :: comp_name
    procedure(set_services_interface) :: set_services
    type(ESMF_Clock), intent(inout) :: clock
    character(*), intent(in) :: namelist_path

    type(comp_info_wrapper_type) info_wrapper
    integer rc

    this%name = comp_name
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

    call ESMF_GridCompInitialize(this%comp, importState=this%import_state, exportState=this%export_state, clock=clock, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    this%import_state = ESMF_StateCreate(name=comp_name, stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    this%export_state = ESMF_StateCreate(name=comp_name, stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

  end subroutine comp_init

  subroutine comp_clear(this)

    class(comp_type), intent(inout) :: this

    integer rc

    call ESMF_GridCompDestroy(this%comp, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Component ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(this%import_state, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Import state ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(this%export_state, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if
    call ESMF_LogWrite('Export state ' // trim(this%name) // ' is destroyed.', ESMF_LOGMSG_INFO)

  end subroutine comp_clear

  subroutine comp_final(this)

    type(comp_type), intent(inout) :: this

    call this%clear()

  end subroutine comp_final

end module comp_wrapper_mod
