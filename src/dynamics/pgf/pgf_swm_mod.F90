! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module pgf_swm_mod

  use const_mod
  use block_mod

  implicit none

contains

  subroutine pgf_swm_prepare(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

  end subroutine pgf_swm_prepare

  subroutine pgf_swm_run(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    integer i, j, k

    associate (mesh => block%mesh, &
               gz   => dstate%gz , & ! in
               du   => dtend%du  , & ! out
               dv   => dtend%dv  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          du%d(i,j,k) = du%d(i,j,k) - (gz%d(i+1,j,k) - gz%d(i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          dv%d(i,j,k) = dv%d(i,j,k) - (gz%d(i,j+1,k) - gz%d(i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine pgf_swm_run

end module pgf_swm_mod
