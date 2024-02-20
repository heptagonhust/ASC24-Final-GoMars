program gmcore_esmf_driver

  use mpi
  use esmf
  use flogger
  use atm_comp_mod
  use lnd_comp_mod
  use perf_mod

  implicit none

  type(comp_type) atm, lnd
  type(ESMF_Clock) clock
  integer thread, rc
  character(256) namelist_path

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call ESMF_InitializePreMPI(rc=rc)

  call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, thread, rc)
  call perf_init()

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, mpiCommunicator=MPI_COMM_WORLD, rc=rc)
  if (rc /= ESMF_SUCCESS) then
    call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  end if
  call ESMF_LogWrite("GMCORE ESMF driver start", ESMF_LOGMSG_INFO)
  call ESMF_LogSet(flush=.true.)

  call atm%init('atm', atm_comp_set_services, clock, namelist_path)
  call lnd%init('lnd', lnd_comp_set_services, clock, namelist_path)

  call atm%run()
  call lnd%run()

  call perf_final()
  call ESMF_Finalize()

end program gmcore_esmf_driver
