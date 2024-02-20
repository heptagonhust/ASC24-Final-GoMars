! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module pgf_mod

  use flogger
  use process_mod, only: proc
  use pgf_swm_mod
  use pgf_lin97_mod
  use pgf_ptb_mod

  implicit none

  private

  public pgf_init
  public pgf_init_after_ic
  public pgf_final
  public pgf_prepare
  public pgf_run

  interface
    subroutine pgf_prepare_interface(block, dstate)
      import block_type, dstate_type
      type(block_type), intent(inout) :: block
      type(dstate_type), intent(inout) :: dstate
    end subroutine pgf_prepare_interface

    subroutine pgf_run_interface(block, dstate, dtend)
      import block_type, dstate_type, dtend_type
      type(block_type), intent(inout) :: block
      type(dstate_type), intent(in) :: dstate
      type(dtend_type), intent(inout) :: dtend
    end subroutine pgf_run_interface
  end interface

  procedure(pgf_prepare_interface), pointer :: pgf_prepare
  procedure(pgf_run_interface), pointer :: pgf_run

contains

  subroutine pgf_init()

    if (baroclinic) then
      select case (pgf_scheme)
      case ('lin97')
        pgf_prepare => pgf_lin97_prepare
        pgf_run => pgf_lin97_run
      case ('ptb')
        pgf_prepare => pgf_ptb_prepare
        pgf_run => pgf_ptb_run
      case default
        if (proc%is_root()) call log_error('Unknown PGF scheme ' // trim(pgf_scheme) // '!')
      end select
    else
      pgf_prepare => pgf_swm_prepare
      pgf_run => pgf_swm_run
    end if

  end subroutine pgf_init

  subroutine pgf_init_after_ic()

    if (baroclinic) then
      select case (pgf_scheme)
      case ('ptb')
        call pgf_ptb_init_after_ic()
      end select
    end if

  end subroutine pgf_init_after_ic

  subroutine pgf_final()

    if (baroclinic) then
      select case (pgf_scheme)
      case ('ptb')
        call pgf_ptb_final()
      end select
    end if

  end subroutine pgf_final

end module pgf_mod
