! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_sfc_mod

  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod

  implicit none

  private

  public mars_nasa_sfc_run

contains

  subroutine mars_nasa_sfc_run(state)

    type(mars_nasa_state_type), intent(inout) :: state

    real(r8) dpt, rib, fh, fm, cdh, cdm, lnz
    integer icol, nlev

    associate (mesh    => state%mesh   , &
               pt      => state%pt     , & ! in
               z       => state%z      , & ! in
               t_sfc   => state%t_sfc  , & ! in
               wsp_bot => state%wsp_bot, & ! in
               z0      => state%z0     , & ! in
               ustar   => state%ustar  , & ! out
               tstar   => state%tstar  )   ! out
    nlev = mesh%nlev
    do icol = 1, mesh%ncol
      dpt = pt(icol,nlev) - t_sfc(icol)
      ! FIXME: Check if z(icol,nlev) is about 5 m.
      ! Calculate bulk Richardson number.
      rib = g * z(icol,nlev) * dpt / (pt(icol,nlev) * wsp_bot(icol)**2 + 1.0e-9_r8)
      ! Calculate stability functions (m for momentum, h for heat).
      if (rib >= 0) then
        fh = 1.0_r8 / (1 + (15 * rib / sqrt(1 + 5 * rib)))
        fm = 1.0_r8 / (1 + (10 * rib / sqrt(1 + 5 * rib)))
      else
        fh = sqrt(1 - 64 * rib)
        fm = sqrt(1 - 16 * rib)
      end if
      lnz = log(z(icol,nlev) / z0(icol))
      ! Calculate drag coefficients.
      cdm = fm * (ka / lnz)**2
      cdh = sqrt(fh) * ka / lnz
      ustar = sqrt(cdm) * wsp_bot(icol)
      tstar = cdh * dpt
    end do
    end associate

  end subroutine mars_nasa_sfc_run

end module mars_nasa_sfc_mod
