! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module pgf_lin97_mod

  use flogger
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use tracer_mod

  implicit none

contains

  subroutine pgf_lin97_prepare(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

  end subroutine pgf_lin97_prepare

  subroutine pgf_lin97_run(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    real(r8) dpk1, dpk2, dgz1, dgz2, dpp1, dpp2, dph1, dph2, L, tmp
    integer i, j, k

    !                    o
    !                   /|
    !                  / |
    !                 /  |
    !                /   |
    !   o-----------/------------o
    !   |          /|            |
    !   |         / |            |
    !   |        /  |            |
    !   |       /   |            |
    !   |      o    |            |
    !   o------|    -------------o
    !          |   /
    !          |  /
    !          | /
    !          |/
    !          o
    associate (mesh    => block%mesh          , & ! in
               qm      => tracers(block%id)%qm, & ! in
               pkh_lev => block%aux%pkh_lev   , & ! in
               ph_lev  => dstate%ph_lev       , & ! in
               gz_lev  => dstate%gz_lev       , & ! in
               p_lev   => dstate%p_lev        , & ! in
               du      => dtend%du            , & ! out
               dv      => dtend%dv            )   ! out
    if (hydrostatic) then
      do k = mesh%full_kds, mesh%full_kde
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i+1,j,k))
            dpk1 = pkh_lev%d(i+1,j,k+1) - pkh_lev%d(i  ,j,k  ) ! 2 - 4
            dpk2 = pkh_lev%d(i  ,j,k+1) - pkh_lev%d(i+1,j,k  ) ! 1 - 3
            dgz1 = gz_lev %d(i  ,j,k+1) - gz_lev %d(i+1,j,k  ) ! 1 - 3
            dgz2 = gz_lev %d(i  ,j,k  ) - gz_lev %d(i+1,j,k+1) ! 4 - 2
            tmp = (dpk1 * dgz1 + dpk2 * dgz2) / mesh%de_lon(j) / (dpk1 + dpk2) / L
            du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dudt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i,j+1,k))
            dpk1 = pkh_lev%d(i,j+1,k+1) - pkh_lev%d(i,j  ,k  ) ! 2 - 4
            dpk2 = pkh_lev%d(i,j  ,k+1) - pkh_lev%d(i,j+1,k  ) ! 1 - 3
            dgz1 = gz_lev %d(i,j  ,k+1) - gz_lev %d(i,j+1,k  ) ! 1 - 3
            dgz2 = gz_lev %d(i,j  ,k  ) - gz_lev %d(i,j+1,k+1) ! 4 - 2
            tmp = (dpk1 * dgz1 + dpk2 * dgz2) / mesh%de_lat(j) / (dpk1 + dpk2) / L
            dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
    else if (nonhydrostatic) then
      do k = mesh%full_kds, mesh%full_kde
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i+1,j,k))
            dpp1 = p_lev  %d(i+1,j,k+1) - p_lev  %d(i  ,j,k  ) - ( &
                   ph_lev %d(i+1,j,k+1) - ph_lev %d(i  ,j,k  )) ! 2 - 4
            dpp2 = p_lev  %d(i  ,j,k+1) - p_lev  %d(i+1,j,k  ) - ( &
                   ph_lev %d(i  ,j,k+1) - ph_lev %d(i+1,j,k  )) ! 1 - 3
            dpk1 = pkh_lev%d(i+1,j,k+1) - pkh_lev%d(i  ,j,k  )  ! 2 - 4
            dpk2 = pkh_lev%d(i  ,j,k+1) - pkh_lev%d(i+1,j,k  )  ! 1 - 3
            dph1 = ph_lev %d(i+1,j,k+1) - ph_lev %d(i  ,j,k  )  ! 2 - 4
            dph2 = ph_lev %d(i  ,j,k+1) - ph_lev %d(i+1,j,k  )  ! 1 - 3
            dgz1 = gz_lev %d(i  ,j,k+1) - gz_lev %d(i+1,j,k  )  ! 1 - 3
            dgz2 = gz_lev %d(i  ,j,k  ) - gz_lev %d(i+1,j,k+1)  ! 4 - 2
            tmp = (                                         &
              (dpk1 * dgz1 + dpk2 * dgz2) / (dpk1 + dpk2) + &
              (dpp1 * dgz1 + dpp2 * dgz2) / (dph1 + dph2)   & ! Nonhydrostatic part
            ) / mesh%de_lon(j) / L
            du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dudt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             |
        !   |             |
        !   |    i,j,k    |
        !   |             |
        !   |             |
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i,j+1,k))
            dpp1 = p_lev  %d(i,j+1,k+1) - p_lev  %d(i,j  ,k  ) - ( &
                   ph_lev %d(i,j+1,k+1) - ph_lev %d(i,j  ,k  )) ! 2 - 4
            dpp2 = p_lev  %d(i,j  ,k+1) - p_lev  %d(i,j+1,k  ) - ( &
                   ph_lev %d(i,j  ,k+1) - ph_lev %d(i,j+1,k  )) ! 1 - 3
            dpk1 = pkh_lev%d(i,j+1,k+1) - pkh_lev%d(i,j  ,k  )  ! 2 - 4
            dpk2 = pkh_lev%d(i,j  ,k+1) - pkh_lev%d(i,j+1,k  )  ! 1 - 3
            dph1 = ph_lev %d(i,j+1,k+1) - ph_lev %d(i,j  ,k  )  ! 2 - 4
            dph2 = ph_lev %d(i,j  ,k+1) - ph_lev %d(i,j+1,k  )  ! 1 - 3
            dgz1 = gz_lev %d(i,j  ,k+1) - gz_lev %d(i,j+1,k  )  ! 1 - 3
            dgz2 = gz_lev %d(i,j  ,k  ) - gz_lev %d(i,j+1,k+1)  ! 4 - 2
            tmp = (                                         &
              (dpk1 * dgz1 + dpk2 * dgz2) / (dpk1 + dpk2) + &
              (dpp1 * dgz1 + dpp2 * dgz2) / (dph1 + dph2)   & ! Nonhydrostatic part
            ) / mesh%de_lat(j) / L
            dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
    end if
    end associate

  end subroutine pgf_lin97_run

end module pgf_lin97_mod
