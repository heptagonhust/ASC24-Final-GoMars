module test_forcing_mod

  use namelist_mod
  use vortex_erosion_test_mod
  use held_suarez_test_mod
  use block_mod

  private

  public test_forcing_run

contains

  subroutine test_forcing_run(dt, time_idx)

    real(r8), intent(in) :: dt
    integer, intent(in) :: time_idx

    integer iblk

    do iblk = 1, size(blocks)
      select case (test_case)
      case ('vortex_erosion')
        call vortex_erosion_test_apply_forcing(blocks(iblk), blocks(iblk)%static)
      case ('held_suarez')
        call held_suarez_test_apply_forcing(blocks(iblk), dt, blocks(iblk)%dstate(time_idx))
      end select
    end do

  end subroutine test_forcing_run

end module test_forcing_mod
