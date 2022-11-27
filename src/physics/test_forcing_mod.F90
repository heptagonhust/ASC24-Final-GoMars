module test_forcing_mod

  use namelist_mod
  use vortex_erosion_test_mod
  use held_suarez_test_mod
  use block_mod

  private

  public test_forcing_run

contains

  subroutine test_forcing_run(block, dt, static, dstate)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(static_type), intent(inout) :: static
    type(dstate_type), intent(inout) :: dstate

    select case (test_case)
    case ('vortex_erosion')
      call vortex_erosion_test_apply_forcing(block, static)
    case ('held_suarez')
      call held_suarez_test_apply_forcing(block, dt, dstate)
    end select

  end subroutine test_forcing_run

end module test_forcing_mod
