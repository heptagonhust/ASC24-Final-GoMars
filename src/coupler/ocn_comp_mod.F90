module ocn_comp_mod

  use mpi
  use mct_mod

  implicit none

  private

  public ocn_comp_init
  public ocn_comp_final

  integer :: comp_id = -1
  integer :: comp_comm = MPI_COMM_NULL

contains

  subroutine ocn_comp_init(comp_id_in, comp_comm_in)

    integer, intent(in) :: comp_id_in
    integer, intent(in) :: comp_comm_in

    comp_id = comp_id_in
    comp_comm = comp_comm_in

  end subroutine ocn_comp_init

  subroutine ocn_comp_final()

  end subroutine ocn_comp_final

end module ocn_comp_mod