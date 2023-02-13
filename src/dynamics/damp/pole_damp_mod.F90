module pole_damp_mod

  use namelist_mod
  use time_mod
  use block_mod
  use filter_mod
  use operators_mod
  use parallel_mod

  implicit none

  private

  public pole_damp_run

contains

  subroutine pole_damp_run(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    if (baroclinic) then
      do k = block%mesh%full_kds, block%mesh%full_kde
        do j = block%mesh%full_jds, block%mesh%full_jde
          do i = block%mesh%full_ids, block%mesh%full_ide
            dstate%pt(i,j,k) = dstate%pt(i,j,k) * dstate%m(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true., &
                     south_halo=.false., north_halo=.false.)
      call filter_on_cell(block%small_filter, dstate%pt)
      do k = block%mesh%full_kds, block%mesh%full_kde
        do j = block%mesh%full_jds, block%mesh%full_jde
          do i = block%mesh%full_ids, block%mesh%full_ide
            dstate%pt(i,j,k) = dstate%pt(i,j,k) / dstate%m(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%filter_halo, dstate%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine pole_damp_run

end module pole_damp_mod
