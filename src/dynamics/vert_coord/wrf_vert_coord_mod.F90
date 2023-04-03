! This vertical coordinate calculation is copied from WRF v4.4.1.
!
! History:
! 
!   - 2021-04-02: First version.

module wrf_vert_coord_mod

  use flogger
  use const_mod
  use namelist_mod, only: ptop

  implicit none

  private

  public wrf_compute_eta

contains

  subroutine wrf_compute_eta(nlev, eta_lev)

    integer, intent(in) :: nlev
    real(r8), intent(out) :: eta_lev(nlev+1)

    real(r8), parameter :: t0           = 288.0_r8    ! Isothermal temperature (K)
    real(r8), parameter :: dzmax        = 3000.0_r8   ! Maximum layer thickness (m)
    real(r8), parameter :: dzbot        = 50.0_r8     ! Bottom layer thickness (m)
    real(r8), parameter :: dzstretch_s  = 1.3_r8      ! Stretching factor for surface layer
    real(r8), parameter :: dzstretch_u  = 1.1_r8      ! Stretching factor for upper layers

    real(r8) ztop, pbot, dz
    real(r8) z_lev(nlev+1), p_lev(nlev+1)
    integer k, k0

    ztop = rd * t0 / g * log(p0 / ptop)
    pbot = p0 * exp(-g * dzbot / (rd * t0))

    eta_lev(1) = 1.0_r8
    eta_lev(2) = (pbot - ptop) / (p0 - ptop)
    z_lev  (1) = 0.0_r8
    z_lev  (2) = dzbot
    p_lev  (1) = p0
    p_lev  (2) = pbot

    ! The stretching transitions from dzstretch_s to dzstretch_u by the time the thickness reaches dzmax/2.
    dz = dzbot
    k0 = 0
    do k = 3, nlev + 1
      dz = dz * dzstretch_u + (dzstretch_s - dzstretch_u) * max((dzmax / 2 - dz) / (dzmax / 2), 0.0_r8)
      if ((ztop - z_lev(k-1)) / (nlev - k + 1) < dz) then
        k0 = k - 1
        exit
      end if
      z_lev(k) = z_lev(k-1) + dz
      p_lev(k) = p0 * exp(-g * z_lev(k) / (rd * t0))
      eta_lev(k) = (p_lev(k) - ptop) / (p0 - ptop)
      if (k == nlev + 1) then
        call log_error('Too few vertical levels for given parameters. Increase nlev!', __FILE__, __LINE__)
      end if
    end do

    dz = (ztop - z_lev(k0)) / (nlev - k0 + 1)
    if (dz > 1.5 * dzmax) then
      call log_error('Upper levels may be too coarse!', __FILE__, __LINE__)
    end if

    do k = k0 + 1, nlev + 1
      z_lev(k) = z_lev(k-1) + dz
      p_lev(k) = p0 * exp(-g * z_lev(k) / (rd * t0))
      eta_lev(k) = (p_lev(k) - ptop) / (p0 - ptop)
    end do
    eta_lev(nlev+1) = 0.0_r8

    ! Reverse order to up-bottom.
    eta_lev = eta_lev(nlev+1:1:-1)

  end subroutine wrf_compute_eta

end module wrf_vert_coord_mod
