module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use filter_mod
  use div_damp_mod
  use vor_damp_mod
  use smag_damp_mod
  use pole_damp_mod
  use laplace_damp_mod

  implicit none

  private

  public damp_init
  public damp_final
  public damp_run

contains

  subroutine damp_init()

    call div_damp_init()
    call vor_damp_init()
    call smag_damp_init()
    call pole_damp_init()
    call laplace_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call div_damp_final()
    call vor_damp_final()
    call smag_damp_final()
    call pole_damp_final()
    call laplace_damp_final()

  end subroutine damp_final

  subroutine damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer cyc

    if (use_pole_damp) then
      call pole_damp_run(block, dstate, dt)
    end if
    if (use_vor_damp) then
      do cyc = 1, vor_damp_cycles
        call vor_damp_run(block, dstate, dt / vor_damp_cycles)
      end do
    end if
    if (use_div_damp) then
      do cyc = 1, div_damp_cycles
        call div_damp_run(block, dstate, dt / div_damp_cycles)
      end do
    end if
    if (use_smag_damp) then
      do cyc = 1, smag_damp_cycles
        call smag_damp_run(block, dstate, dt / smag_damp_cycles)
      end do
    end if

  end subroutine damp_run

end module damp_mod
