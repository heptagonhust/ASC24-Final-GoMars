module drydep_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: gravit, rair
  use ppgrid

  implicit none

contains

  subroutine inidrydep(xrair, xgravit)

    ! Initialize dry deposition parameterization.

    real(r8), intent(in) :: xrair                ! Gas constant for dry air
    real(r8), intent(in) :: xgravit              ! Gravitational acceleration

  end subroutine inidrydep

  subroutine setdvel(ncol, landfrac, icefrac, ocnfrac, vgl, vgo, vgsi, vg)

    ! Set the deposition velocity depending on whether we are over
    ! land, ocean, and snow/ice

    integer, intent(in) :: ncol
    real (r8), intent(in) :: landfrac(pcols)    ! land fraction
    real (r8), intent(in) :: icefrac(pcols)     ! ice fraction
    real (r8), intent(in) :: ocnfrac(pcols)     ! ocean fraction

    real(r8), intent(in) :: vgl                 ! dry deposition velocity in m/s (land)
    real(r8), intent(in) :: vgo                 ! dry deposition velocity in m/s (ocean)
    real(r8), intent(in) :: vgsi                ! dry deposition velocity in m/s (snow/ice)

    ! Output arguments:
    real(r8), intent(out) :: vg(pcols)          ! dry deposition velocity in m/s

    ! Local variables:
    integer i
    real(r8) a

    do i = 1, ncol
      vg(i) = landfrac(i)*vgl + ocnfrac(i)*vgo + icefrac(i)*vgsi
      ! if (ioro(i).eq.0) then
      !   vg(i) = vgo
      ! else if (ioro(i).eq.1) then
      !   vg(i) = vgl
      ! else
      !   vg(i) = vgsi
      ! end if
    end do

  end subroutine setdvel

  subroutine ddflux(ncol, vg, q, p, tv, flux)

    ! Compute surface flux due to dry deposition processes.

    ! Input arguments:
    integer , intent(in) :: ncol
    real(r8), intent(in) :: vg(pcols)  ! dry deposition velocity in m/s
    real(r8), intent(in) :: q(pcols)   ! tracer conc. in surface layer (kg tracer/kg moist air)
    real(r8), intent(in) :: p(pcols)   ! midpoint pressure in surface layer (Pa)
    real(r8), intent(in) :: tv(pcols)  ! midpoint virtual temperature in surface layer (K)

    ! Output arguments:
    real(r8), intent(out) :: flux(pcols) ! flux due to dry deposition in kg/m^s/sec

    ! Local variables:
    integer i

    do i = 1, ncol
      flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
    end do

  end subroutine ddflux

  subroutine d3ddflux(ncol, vlc_dry, q, pmid, pdel, tv, dep_dry, dep_dry_tend, dt)

    ! Do 3d-settling deposition calculations following Zender's dust codes, Dec 02.
    ! Author: Natalie Mahowald

    ! Input arguments:
    integer , intent(in) :: ncol
    real(r8), intent(in) :: vlc_dry(pcols,pver) ! dry deposition velocity in m/s
    real(r8), intent(in) :: q(pcols,pver)       ! tracer conc. in surface layer (kg tracer/kg moist air)
    real(r8), intent(in) :: pmid(pcols,pver)    ! midpoint pressure in surface layer (Pa)
    real(r8), intent(in) :: pdel(pcols,pver)    ! delta pressure across level (Pa)
    real(r8), intent(in) :: tv(pcols,pver)      ! midpoint virtual temperature in surface layer (K)
    real(r8), intent(in) :: dt

    ! Output arguments:
    real(r8), intent(out) :: dep_dry(pcols)           ! flux due to dry deposition in kg /m^s/sec
    real(r8), intent(out) :: dep_dry_tend(pcols,pver) ! flux due to dry deposition in kg /m^s/sec

    ! Local variables:
    real(r8) flux(pcols,0:pver)  ! downward flux at each level:  kg/m2/s
    integer i, k

    do i = 1, ncol
      flux(i,0) = 0._r8
    end do
    do k = 1, pver
      do i = 1, ncol
        flux(i,k) = -min(vlc_dry(i,k) * q(i,k) * pmid(i,k) /(tv(i,k) * rair), &
                      q(i,k)*pdel(i,k)/gravit/dt)
        dep_dry_tend(i,k)=(flux(i,k)-flux(i,k-1))/pdel(i,k)*gravit  !kg/kg/s
      end do
    end do
    ! surface flux:
    do i = 1, ncol
      dep_dry(i) = flux(i,pver)
    end do

  end subroutine d3ddflux

  subroutine  calcram(ncol, landfrac, icefrac, ocnfrac, obklen, &
                      ustar, ram1in, ram1, t, pmid, pdel, fvin, fv)

    ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
    ! from Seinfeld and Pandis, p.963.
    ! Author: Natalie Mahowald

    integer,  intent(in) :: ncol
    real(r8), intent(in) :: ram1in(pcols)
    real(r8), intent(in) :: fvin(pcols)
    real(r8), intent(out) :: ram1(pcols)
    real(r8), intent(out) :: fv(pcols)
    real(r8), intent(in) :: obklen(pcols)
    real(r8), intent(in) :: ustar(pcols)
    real(r8), intent(in) :: landfrac(pcols)
    real(r8), intent(in) :: icefrac(pcols)
    real(r8), intent(in) :: ocnfrac(pcols)
    real(r8), intent(in) :: t(pcols)
    real(r8), intent(in) :: pmid(pcols)
    real(r8), intent(in) :: pdel(pcols)

    real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
    real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
    real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

    real(r8) z, psi, psi0, nu, nu0, temp, ram
    integer i

    do i = 1, ncol
      z = pdel(i) * rair * t(i) / pmid(i) / gravit / 2.0_r8   ! use half the layer height like Ganzefeld and Lelieveld, 1995
      if (obklen(i) == 0) then
        psi  = 0
        psi0 = 0
      else
        psi  = min(max(z / obklen(i), -1.0_r8), 1.0_r8)
        psi0 = min(max(zzocen / obklen(i), -1.0_r8), 1.0_r8)
      end if
      temp = z / zzocen
      if (icefrac(i) > 0.5_r8) then
        if (obklen(i) > 0) then
          psi0 = min(max(zzsice / obklen(i), -1.0_r8), 1.0_r8)
        else
          psi0 = 0
        end if
        temp = z / zzsice
      end if
      if (psi > 0._r8) then
        ram = 1 / xkar / ustar(i) * (log(temp) + 4.7_r8 * (psi - psi0))
      else
        nu  = (1.000_r8 - 15.000_r8 * psi)**(.25_r8)
        nu0 = (1.000_r8 - 15.000_r8 *psi0)**(.25_r8)
        if (ustar(i) /= 0) then
          ram = 1 / xkar / ustar(i) * (log(temp) &
              + log(((nu0**2 + 1.00_r8) * (nu0+1.0_r8)**2) / ((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
              + 2.0_r8*(atan(nu)-atan(nu0)))
        else
          ram = 0
        end if
      end if
      if (landfrac(i) < 0.000000001_r8) then
        fv  (i) = ustar(i)
        ram1(i) = ram
      else
        fv  (i) = fvin(i)
        ram1(i) = ram1in(i)
      end if
    end do

    ! fvitt -- fv == 0 causes a floating point exception in
    ! dry dep of sea salts and dust
    where (fv(:ncol) == 0)
      fv(:ncol) = 1.e-12_r8
    end where

  end subroutine calcram

end module drydep_mod
