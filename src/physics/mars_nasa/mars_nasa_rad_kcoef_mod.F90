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
! Introduction:
!
!   The absorption coefficient K is a function of local temperature, pressure, 
!   gas mixture, and wavelength. It is pre-computed over a wide range of air
!   conditions, and stored in a look up table.
!
!   K distribution method (KDM) resorts K in each selected spectral band by
!   assuming each wavelength is no more important than others.
!
! ==============================================================================

module mars_nasa_rad_kcoef_mod

  use fiona
  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_spectra_mod

  implicit none

  private

  public mars_nasa_rad_kcoef_init
  public mars_nasa_rad_kcoef_final
  public get_kcoef_vis
  public ntref
  public npint
  public nqref
  public ngauss
  public tref
  public pref
  public logpint
  public qh2oref
  public qco2ref
  public gwgt
  public klut_vis
  public klut_ir
  public f0_vis
  public f0_ir

  ! Number of reference temperature levels for K-coefficients
  integer ntref
  ! Number of reference pressure levels for K-coefficients
  integer npref
  ! Number of Lagrange interpolated reference pressures for the CO2 K-coefficients
  integer npint
  ! Number of different water mixing ratios for K-coefficients that are now CO2 + H2O
  integer nqref
  ! Number of Gaussian quadrature points for K-coefficients
  integer ngauss

  ! Reference temperature for K-coefficients (K)
  real(r8), allocatable, dimension(:        ) :: tref
  ! Reference pressure for K-coefficients (Pa)
  real(r8), allocatable, dimension(  :      ) :: pref
  ! Interpolated reference pressure for K-coefficients (logPa)
  real(r8), allocatable, dimension(  :      ) :: logpint
  ! Reference water vapor volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(    :    ) :: qh2oref
  ! Reference CO2 volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(    :    ) :: qco2ref
  ! Gaussian point weights
  real(r8), allocatable, dimension(        :) :: gwgt
  ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_vis
  ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_ir
  ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_vis
  ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_ir
  ! Fraction of zeros in visible CO2 K-coefficients
  real(r8), allocatable, dimension(      :  ) :: f0_vis
  ! Fraction of zeros in infrared C02 K-coefficients
  real(r8), allocatable, dimension(      :  ) :: f0_ir

contains

  subroutine mars_nasa_rad_kcoef_init()

    call mars_nasa_rad_kcoef_final()

    call fiona_open_dataset('kcoef', file_path=kcoef_file)
    call fiona_get_dim('kcoef', 'tref' , size=ntref )
    call fiona_get_dim('kcoef', 'pref' , size=npref )
    call fiona_get_dim('kcoef', 'qref' , size=nqref )
    call fiona_get_dim('kcoef', 'gauss', size=ngauss)
    allocate(tref       (ntref                              ))
    allocate(pref       (      npref                        ))
    allocate(qco2ref    (            nqref                  ))
    allocate(qh2oref    (            nqref                  ))
    allocate(gwgt       (                             ngauss))
    allocate(klut_in_vis(ntref,npref,nqref,spec_vis%n,ngauss))
    allocate(klut_in_ir (ntref,npref,nqref,spec_ir %n,ngauss))
    allocate(f0_vis     (                  spec_vis%n       ))
    allocate(f0_ir      (                  spec_ir %n       ))
    call fiona_start_input('kcoef')
    call fiona_input('kcoef', 'tref'    , tref       )
    call fiona_input('kcoef', 'pref'    , pref       )
    call fiona_input('kcoef', 'qco2ref' , qco2ref    )
    call fiona_input('kcoef', 'qh2oref' , qh2oref    )
    call fiona_input('kcoef', 'gwgt'    , gwgt       )
    call fiona_input('kcoef', 'klut_vis', klut_in_vis)
    call fiona_input('kcoef', 'klut_ir' , klut_in_ir )
    call fiona_input('kcoef', 'f0_vis'  , f0_vis     )
    call fiona_input('kcoef', 'f0_ir'   , f0_ir      )
    call fiona_end_input('kcoef')

    call interp_kcoef()

  end subroutine mars_nasa_rad_kcoef_init

  subroutine mars_nasa_rad_kcoef_final()

    if (allocated(tref       )) deallocate(tref       )
    if (allocated(pref       )) deallocate(pref       )
    if (allocated(logpint    )) deallocate(logpint    )
    if (allocated(qco2ref    )) deallocate(qco2ref    )
    if (allocated(qh2oref    )) deallocate(qh2oref    )
    if (allocated(gwgt       )) deallocate(gwgt       )
    if (allocated(klut_in_vis)) deallocate(klut_in_vis)
    if (allocated(klut_in_ir )) deallocate(klut_in_ir )
    if (allocated(klut_vis   )) deallocate(klut_vis   )
    if (allocated(klut_ir    )) deallocate(klut_ir    )
    if (allocated(f0_vis     )) deallocate(f0_vis     )
    if (allocated(f0_ir      )) deallocate(f0_ir      )

  end subroutine mars_nasa_rad_kcoef_final

  subroutine interp_kcoef()

    real(r8) logpref(npref)
    real(r8) p1, p2, dp
    integer it, ip, iq, ig, iw, i, ipp

    ! Take log10 of the reference pressures.
    do ip = 1, npref
      logpref(ip) = log10(pref(ip))
    end do

    ! Refine each reference pressure bin into 5 smaller bins.
    npint = (npref - 1) * 5 + 1

    allocate(logpint(npint))
    do ip = 1, npref - 1
      p1 = logpref(ip)
      p2 = logpref(ip+1)
      dp = (p2 - p1) / 5.0_r8
      do i = 1, 5
        ipp = (ip - 1) * 5 + i
        logpint(ipp) = p1 + (i - 1) * dp
      end do
    end do
    logpint(npint) = logpref(npref)

    ! Take log10 of the values, since the smallest value is 1.0e-200.
    do ig = 1, ngauss
    do iq = 1, nqref
    do ip = 1, npref
    do it = 1, ntref
      do iw = 1, spec_vis%n
        if (klut_in_vis(it,ip,iq,iw,ig) > 1.0d-200) then
          klut_in_vis(it,ip,iq,iw,ig) = log10(klut_in_vis(it,ip,iq,iw,ig))
        else
          klut_in_vis(it,ip,iq,iw,ig) = -200
        end if
      end do
      do iw = 1, spec_ir%n
        if (klut_in_ir (it,ip,iq,iw,ig) > 1.0d-200) then
          klut_in_ir (it,ip,iq,iw,ig) = log10(klut_in_ir (it,ip,iq,iw,ig))
        else
          klut_in_ir (it,ip,iq,iw,ig) = -200
        end if
      end do
    end do
    end do
    end do
    end do

    allocate(klut_vis(ntref,npint,nqref,spec_vis%n,ngauss))
    allocate(klut_ir (ntref,npint,nqref,spec_ir%n ,ngauss))

    do ig = 1, ngauss
    do iw = 1, spec_vis%n
    do iq = 1, nqref
    do it = 1, ntref
      ip = 1
      do i = 1, 5
        ipp = (ip - 1) * 5 + i
        klut_vis(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
          logpref(ip:ip+3), klut_in_vis(it,ip:ip+3,iq,iw,ig), logpint(ipp))
      end do
      do ip = 2, npref - 2
        do i = 1, 5
          ipp = (ip - 1) * 5 + i
          klut_vis(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
            logpref(ip-1:ip+2), klut_in_vis(it,ip-1:ip+2,iq,iw,ig), logpint(ipp))
        end do
      end do
      ip = npref - 1
      do i = 1, 5
        ipp = (ip - 1) * 5 + i
        klut_vis(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
          logpref(ip-2:ip+1), klut_in_vis(it,ip-2:ip+1,iq,iw,ig), logpint(ipp))
      end do
      klut_vis(it,npint,iq,iw,ig) = 10 ** klut_in_vis(it,npref,iq,iw,ig)
    end do
    end do
    end do
    end do

    do ig = 1, ngauss
    do iw = 1, spec_ir%n
    do iq = 1, nqref
    do it = 1, ntref
      ip = 1
      do i = 1, 5
        ipp = (ip - 1) * 5 + i
        klut_ir(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
          logpref(ip:ip+3), klut_in_ir(it,ip:ip+3,iq,iw,ig), logpint(ipp))
      end do
      do ip = 2, npref - 2
        do i = 1, 5
          ipp = (ip - 1) * 5 + i
          klut_ir(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
            logpref(ip-1:ip+2), klut_in_ir(it,ip-1:ip+2,iq,iw,ig), logpint(ipp))
        end do
      end do
      ip = npref - 1
      do i = 1, 5
        ipp = (ip - 1) * 5 + i
        klut_ir(it,ipp,iq,iw,ig) = 10 ** lagrange_interp_4( &
          logpref(ip-2:ip+1), klut_in_ir(it,ip-2:ip+1,iq,iw,ig), logpint(ipp))
      end do
      klut_ir(it,npint,iq,iw,ig) = 10 ** klut_in_ir(it,npref,iq,iw,ig)
    end do
    end do
    end do
    end do

  contains

    pure real(r8) function lagrange_interp_4(xi, yi, x) result(y)

      real(r8), intent(in) :: xi(4)
      real(r8), intent(in) :: yi(4)
      real(r8), intent(in) :: x

      real(r8) x1, x2, x3, x4

      x1 = x - xi(1)
      x2 = x - xi(2)
      x3 = x - xi(3)
      x4 = x - xi(4)

      y = x2 * x3 * x4 * yi(1) / ((xi(1) - xi(2)) * (xi(1) - xi(3)) * (xi(1) - xi(4))) + &
          x1 * x3 * x4 * yi(2) / ((xi(2) - xi(1)) * (xi(2) - xi(3)) * (xi(2) - xi(4))) + &
          x1 * x2 * x4 * yi(3) / ((xi(3) - xi(1)) * (xi(3) - xi(2)) * (xi(3) - xi(4))) + &
          x1 * x2 * x3 * yi(4) / ((xi(4) - xi(1)) * (xi(4) - xi(2)) * (xi(4) - xi(3)))

    end function lagrange_interp_4

  end subroutine interp_kcoef

  subroutine get_kcoef_vis(t, p, qh2o, kcoef)

    ! Temperature (K)
    real(r8), intent(in   ) :: t
    ! Pressure (Pa)
    real(r8), intent(in   ) :: p
    ! Water vapor mixing ratio (kg kg-1)
    real(r8), intent(in   ) :: qh2o
    ! Optical depth
    real(r8), intent(inout) :: kcoef(:,:)

    integer i, it, ip, iq, is, ig
    real(r8) logp, wt, wp, wq, c(4), k(4)

    if (t < tref(1)) then
      it = 1
    else if (t > tref(ntref)) then
      it = ntref - 1
    else
      do i = 1, ntref - 1
        if (tref(i) < t .and. t <= tref(i+1)) then
          it = i
          exit
        end if
      end do
    end if
    wt = (t - tref(it)) / (tref(it+1) - tref(it))

    logp = log10(p)
    if (logp < logpint(1)) then
      ip = 1
    else if (logp > logpint(npint)) then
      ip = npint - 1
    else
      do i = 1, npint - 1
        if (logpint(i) < logp .and. logp <= logpint(i+1)) then
          ip = i
          exit
        end if
      end do
    end if
    wp = (logp - logpint(ip)) / (logpint(ip+1) - logpint(ip))

    c(1) = (1 - wp) * (1 - wt)
    c(2) = wp * (1 - wt)
    c(3) = wp * wt
    c(4) = (1 - wp) * wt

    if (qh2o <= qh2oref(1)) then
      iq = 1
      wq = 0
    else if (qh2o >= qh2oref(nqref)) then
      iq = nqref
      wq = 0
    else
      do i = 2, nqref
        if (qh2oref(i-1) <= qh2o .and. qh2o < qh2oref(i)) then
          iq = i - 1
          wq = (qh2o - qh2oref(i-1)) / (qh2oref(i) - qh2oref(i-1))
          exit
        end if
      end do
    end if

    do is = 1, spec_vis%n
      do ig = 1, ngauss - 1
        k(1) = klut_vis(it  ,ip  ,iq,is,ig) + wq * (klut_vis(it  ,ip  ,iq+1,is,ig) - klut_vis(it  ,ip  ,iq,is,ig))
        k(2) = klut_vis(it  ,ip+1,iq,is,ig) + wq * (klut_vis(it  ,ip+1,iq+1,is,ig) - klut_vis(it  ,ip+1,iq,is,ig))
        k(3) = klut_vis(it+1,ip+1,iq,is,ig) + wq * (klut_vis(it+1,ip+1,iq+1,is,ig) - klut_vis(it+1,ip+1,iq,is,ig))
        k(4) = klut_vis(it+1,ip  ,iq,is,ig) + wq * (klut_vis(it+1,ip  ,iq+1,is,ig) - klut_vis(it+1,ip  ,iq,is,ig))
        kcoef(is,ig) = sum(c * k)
      end do
    end do

  end subroutine get_kcoef_vis

end module mars_nasa_rad_kcoef_mod
