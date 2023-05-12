module ice_comp_mod

  use mpi
  use mct_mod

  implicit none

  private

  public ice_comp_init
  public ice_comp_final

  integer :: comp_id = -1
  integer :: comp_comm = MPI_COMM_NULL

contains

  subroutine ice_comp_init(comp_id_in, comp_comm_in)

    integer, intent(in) :: comp_id_in
    integer, intent(in) :: comp_comm_in

    comp_id = comp_id_in
    comp_comm = comp_comm_in

  end subroutine ice_comp_init

  subroutine ice_comp_final()

  end subroutine ice_comp_final

end module ice_comp_mod