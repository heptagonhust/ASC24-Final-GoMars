module pgf_swm_mod

  use const_mod
  use parallel_mod
  use block_mod

  implicit none

contains

  subroutine pgf_swm_prepare(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

  end subroutine pgf_swm_prepare

  subroutine pgf_swm_run(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    integer i, j, k

    associate (mesh    => block%mesh   , &
               gz      => dstate%gz    , & ! in
               pgf_lon => dtend%pgf_lon, & ! out
               pgf_lat => dtend%pgf_lat)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pgf_lon(i,j,k) = (gz(i+1,j,k) - gz(i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pgf_lat(i,j,k) = (gz(i,j+1,k) - gz(i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine pgf_swm_run

end module pgf_swm_mod
