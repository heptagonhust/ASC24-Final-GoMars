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
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          pgf_lon(i,j,k) = (gz(i+1,j,k) - gz(i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          pgf_lat(i,j,k) = (gz(i,j+1,k) - gz(i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine pgf_swm_run

end module pgf_swm_mod
