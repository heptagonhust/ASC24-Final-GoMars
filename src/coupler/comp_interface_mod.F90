module comp_interface_mod

  implicit none

  interface
    subroutine comp_init_interface(comp_id, comp_comm)
      integer, intent(in) :: comp_id
      integer, intent(in) :: comp_comm
    end subroutine comp_init_interface
    subroutine comp_final_interface()
    end subroutine comp_final_interface
  end interface

  type comp_agent_type
    procedure(comp_init_interface), pointer, nopass :: init => null()
    procedure(comp_final_interface), pointer, nopass :: final => null()
  end type comp_agent_type

end module comp_interface_mod