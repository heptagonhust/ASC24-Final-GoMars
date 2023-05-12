module atm_comp_mod

  use mpi
  use mct_mod

  implicit none

  private

  public atm_comp_init
  public atm_comp_final

  integer :: comp_id = -1
  integer :: comp_comm = MPI_COMM_NULL

contains

  subroutine atm_comp_init(comp_id_in, comp_comm_in)

    integer, intent(in) :: comp_id_in
    integer, intent(in) :: comp_comm_in

    comp_id = comp_id_in
    comp_comm = comp_comm_in

  end subroutine atm_comp_init

  subroutine atm_comp_final()

  end subroutine atm_comp_final

end module atm_comp_mod