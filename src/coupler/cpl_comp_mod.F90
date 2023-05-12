module cpl_comp_mod

  use mpi
  use mct_mod
  use comp_interface_mod
  use atm_comp_mod
  use lnd_comp_mod
  use ocn_comp_mod
  use ice_comp_mod

  implicit none

  private

  logical :: do_manage_mpi = .false.
  integer :: ncomps = 0
  integer :: global_comm = MPI_COMM_NULL
  integer, pointer :: comp_ids(:)
  integer, pointer :: comp_comms(:)

  type(comp_agent_type), allocatable :: comp_agents(:)

contains

  subroutine cpl_comp_init(global_comm_in, active_atm, active_lnd, active_ocn, active_ice)

    integer, intent(in), optional :: global_comm_in
    logical, intent(in), optional :: active_atm
    logical, intent(in), optional :: active_lnd
    logical, intent(in), optional :: active_ocn
    logical, intent(in), optional :: active_ice

    integer i, ierr

    if (merge(active_atm, .true. , present(active_atm))) ncomps = ncomps + 1
    if (merge(active_lnd, .false., present(active_lnd))) ncomps = ncomps + 1
    if (merge(active_ocn, .false., present(active_ocn))) ncomps = ncomps + 1
    if (merge(active_ice, .false., present(active_ice))) ncomps = ncomps + 1

    global_comm = merge(global_comm_in, MPI_COMM_WORLD, present(global_comm_in))
    if (.not. present(global_comm_in)) then
      call MPI_INIT(ierr)
      do_manage_mpi = .true.
    end if

    allocate(comp_ids  (ncomps))
    allocate(comp_comms(ncomps))
    allocate(comp_agents(ncomps))

    do i = 1, ncomps
      comp_ids(i) = i
      call MPI_COMM_DUP(global_comm, comp_comms(i), ierr)
    end do

    call mct_world_init(ncomps, global_comm, comp_comms, comp_ids)

    i = 1
    if (merge(active_atm, .true. , present(active_atm))) then
      comp_agents(i)%init  => atm_comp_init
      comp_agents(i)%final => atm_comp_final
      i = i + 1
    end if
    if (merge(active_lnd, .false., present(active_lnd))) then
      comp_agents(i)%init  => lnd_comp_init
      comp_agents(i)%final => lnd_comp_final
      i = i + 1
    end if
    if (merge(active_ocn, .false., present(active_ocn))) then
      comp_agents(i)%init  => ocn_comp_init
      comp_agents(i)%final => ocn_comp_final
      i = i + 1
    end if
    if (merge(active_ice, .false., present(active_ice))) then
      comp_agents(i)%init  => ice_comp_init
      comp_agents(i)%final => ice_comp_final
      i = i + 1
    end if

    do i = 1, ncomps
      call comp_agents(i)%init(comp_ids(i), comp_comms(i))
    end do

  end subroutine cpl_comp_init

  subroutine cpl_comp_final()

    integer i, ierr

    do i = 1, ncomps
      call comp_agents(i)%final()
    end do

    if (associated(comp_ids  )) deallocate(comp_ids  )
    if (associated(comp_comms)) deallocate(comp_comms)

    if (do_manage_mpi) then
      call MPI_FINALIZE(ierr)
    end if

  end subroutine cpl_comp_final

end module cpl_comp_mod
