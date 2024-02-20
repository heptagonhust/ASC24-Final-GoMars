! This module wraps GPTL perf library.

module perf_mod

#ifdef HAS_GPTL
  use gptl
#endif
  use mpi

  implicit none

  public perf_init
  public perf_start
  public perf_stop
  public perf_final

  ! Just for CAM
  public t_startf
  public t_stopf
  public t_disablef
  public t_enablef
  public t_adj_detailf
  public t_barrierf

  logical :: initialized    = .false.
  logical :: barrier        = .false.
  integer :: disable_depth  = 0
  integer :: detail_limit   = 1
  integer :: curr_detail    = 0

contains

  subroutine perf_init()

    integer ierr

#ifdef HAS_GPTL
    ierr = gptlsetoption(gptloverhead, 0)
    ierr = gptlsetoption(gptlpercent, 0)
    ierr = gptlsetoption(gptlabort_on_error, 1)
    ierr = gptlsetoption(gptlsync_mpi, 1)
    if (gptlinitialize() < 0) then
      stop 'Failed to initialize GPTL!'
    end if
    ierr = gptlstart('total')
#endif

    initialized = .true.

  end subroutine perf_init

  subroutine perf_start(name)

    character(*), intent(in) :: name

    integer ierr

#ifdef HAS_GPTL
    ierr = gptlstart(name)
#endif

  end subroutine perf_start

  subroutine perf_stop(name)

    character(*), intent(in) :: name

    integer ierr

#ifdef HAS_GPTL
    ierr = gptlstop(name)
#endif

  end subroutine perf_stop

  subroutine perf_final()

    integer ierr

#ifdef HAS_GPTL
    ierr = gptlstop('total')
    ierr = gptlpr(0)
    if (gptlfinalize() /= 0) then
      stop 'Failed to call finalize GPTL!'
    end if
#endif

  end subroutine perf_final

  subroutine t_startf(event)

    character(*), intent(in) :: event

#ifdef HAS_GPTL
    if (GPTLstart(event) /= 0) then
      stop 'Failed to call GPTLstart!'
    end if
#endif

  end subroutine t_startf

  subroutine t_stopf(event)

    character(*), intent(in) :: event

#ifdef HAS_GPTL
    if (GPTLstop(event) /= 0) then
      stop 'Failed to call GPTLstop!'
    end if
#endif

  end subroutine t_stopf

  subroutine t_disablef()

#ifdef HAS_GPTL
    if (.not. initialized) return
    if (disable_depth == 0) then
      if (GPTLdisable() /= 0) then
        stop 'Failed to call GPTLdisable!'
      end if
    end if
    disable_depth = disable_depth + 1
#endif

  end subroutine t_disablef

  subroutine t_enablef()

#ifdef HAS_GPTL
    if (.not. initialized) return
    if (disable_depth > 0) then
      if (disable_depth == 1) then
        if (GPTLenable() /= 0) then
          stop 'Failed to call GPTLenable!'
        end if
      end if
      disable_depth = disable_depth - 1
    end if
#endif

  end subroutine t_enablef

  subroutine t_adj_detailf(detail_adjustment)

    integer, intent(in) :: detail_adjustment

#ifdef HAS_GPTL
    if (.not. initialized) return
    ! Use disable/enable to implement timing_detail logic so also control direct GPTL calls.
    if (curr_detail <= detail_limit .and. curr_detail + detail_adjustment > detail_limit) then
      call t_disablef()
    else if (curr_detail > detail_limit .and. curr_detail + detail_adjustment <= detail_limit) then
      call t_enablef()
    end if

    curr_detail = curr_detail + detail_adjustment
#endif

  end subroutine t_adj_detailf

  subroutine t_barrierf(event, comm)

    character(*), intent(in), optional :: event
    integer, intent(in), optional :: comm

    integer ierr

#ifdef HAS_GPTL
    if (.not. initialized) return
    if (disable_depth > 0) return

    if (barrier) then
      if (present(event)) call t_startf(event)
      if (present(comm)) then
        call MPI_Barrier(comm, ierr)
      else
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end if
      if (present(event)) call t_stopf(event)
    end if
#endif

  end subroutine t_barrierf

end module perf_mod
