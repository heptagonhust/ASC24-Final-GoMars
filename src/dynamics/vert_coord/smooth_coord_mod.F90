module smooth_coord_mod

  use const_mod
  use namelist_mod
  use wrf_vert_coord_mod
  use latlon_mesh_mod, only: global_mesh

  implicit none

  private

  public smooth_coord_init
  public smooth_coord_final
  public smooth_coord_calc_mg
  public smooth_coord_calc_mg_lev
  public smooth_coord_calc_dmgdt
  public smooth_coord_calc_dmgdt_lev
  public smooth_coord_calc_ddmgdt
  public nlevp

  real(r8), parameter :: eta_b = 0.2_r8
  real(r8), parameter :: eta_c = 0.8_r8
  real(r8), allocatable :: bi(:), ci(:)
  real(r8), allocatable :: bm(:), cm(:)
  integer :: nlevp = 0

contains

  subroutine smooth_coord_init()

    real(r8) b0, b1, b2, b3
    real(r8) c0, c1, c2, c3
    real(r8) eta
    integer k

    if (allocated(bi)) deallocate(bi)
    if (allocated(ci)) deallocate(ci)
    if (allocated(bm)) deallocate(bm)
    if (allocated(cm)) deallocate(cm)

    allocate(bi(nlev+1)); bi = 0
    allocate(ci(nlev+1)); ci = 0
    allocate(bm(nlev  )); bm = 0
    allocate(cm(nlev  )); cm = 0

    call wrf_compute_eta(nlev, global_mesh%half_lev)

    b0 = 2 * eta_b**2 / (1 - eta_b)**3
    b1 = -eta_b * (4 + eta_b + eta_b**2) / (1 - eta_b)**3
    b2 = 2 * (1 + eta_b + eta_b**2) / (1 - eta_b)**3
    b3 = -(1 + eta_b) / (1 - eta_b)**3
    c0 = 2 * eta_c**2 / (1 - eta_c)**3
    c1 = -eta_c * (4 + eta_c + eta_c**2) / (1 - eta_c)**3
    c2 = 2 * (1 + eta_c + eta_c**2) / (1 - eta_c)**3
    c3 = -(1 + eta_c) / (1 - eta_c)**3
    do k = 1, nlev + 1
      eta = global_mesh%half_lev(k)
      if (eta >= eta_b) bi(k) = b0 + b1 * eta + b2 * eta**2 + b3 * eta**3
      if (eta >= eta_c) ci(k) = c0 + c1 * eta + c2 * eta**2 + c3 * eta**3
    end do

    do k = 1, nlev
      bm(k) = (bi(k) + bi(k+1)) * 0.5_r8
      cm(k) = (ci(k) + ci(k+1)) * 0.5_r8
    end do

    do k = 1, nlev
      global_mesh%full_lev(k) = (global_mesh%half_lev(k) + global_mesh%half_lev(k+1)) * 0.5_r8
    end do

    nlevp = 0
    do k = 1, nlev
      if (bm(k) == 0 .and. cm(k) == 0) nlevp = nlevp + 1
    end do

  end subroutine smooth_coord_init

  subroutine smooth_coord_final()

    if (allocated(bi)) deallocate(bi)
    if (allocated(ci)) deallocate(ci)
    if (allocated(bm)) deallocate(bm)
    if (allocated(cm)) deallocate(cm)

  end subroutine smooth_coord_final

  pure real(r8) function smooth_coord_calc_mg(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    real(r8) eta

    eta = global_mesh%full_lev(k)

    res = bm(k) * (mgs - ref_ps_perb - ptop) + (eta - bm(k)) * (p0 - ptop) + cm(k) * ref_ps_perb + ptop

  end function smooth_coord_calc_mg

  pure real(r8) function smooth_coord_calc_mg_lev(k, mgs, ref_ps_perb) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: mgs
    real(r8), intent(in), optional :: ref_ps_perb

    real(r8) eta

    eta = global_mesh%half_lev(k)

    res = bi(k) * (mgs - ref_ps_perb - ptop) + (eta - bi(k)) * (p0 - ptop) + ci(k) * ref_ps_perb + ptop

  end function smooth_coord_calc_mg_lev

  pure real(r8) function smooth_coord_calc_dmgdt(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = bm(k) * dmgsdt

  end function smooth_coord_calc_dmgdt

  pure real(r8) function smooth_coord_calc_dmgdt_lev(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = bi(k) * dmgsdt

  end function smooth_coord_calc_dmgdt_lev

  pure real(r8) function smooth_coord_calc_ddmgdt(k, dmgsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dmgsdt

    res = (bi(k+1) - bi(k)) * dmgsdt

  end function smooth_coord_calc_ddmgdt

end module smooth_coord_mod
