module mars_nasa_rad_mod

  use fiona
  use const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_spectra_mod
  use mars_nasa_optics_mod

  implicit none

  integer               nlev_rad
  integer               ntref          ! Number of reference temperature levels for K-coefficients
  integer               npref          ! Number of reference pressure levels for K-coefficients
  integer               nqref          ! Number of different water mixing ratios for K-coefficients that are now CO2 + H2O
  integer               ngauss         ! Number of Gaussian quadrature points for K-coefficients
  integer               npint          ! Number of Lagrange interpolated reference pressures for the CO2 K-coefficients
  integer               nbin_dust      ! Number of dust particle size bins
  integer, parameter :: taumax    = 35 ! Maximum optical depth

  real(r8), parameter :: ray_p0   = 9.423e6_r8        ! Reference pressure for Rayleigh scattering (Pa)

  real(r8), allocatable, dimension(:        ) :: tref        ! Reference temperature for K-coefficients (K)
  real(r8), allocatable, dimension(  :      ) :: pref        ! Reference pressure for K-coefficients (Pa)
  real(r8), allocatable, dimension(    :    ) :: qh2oref     ! Reference water vapor volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(    :    ) :: qco2ref     ! Reference CO2 volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(        :) :: gwgt        ! Gaussian point weights
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_vis ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_ir  ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_vis    ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_ir     ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(      :  ) :: f0_vis      ! Fraction of zeros in visible CO2 K-coefficients
  real(r8), allocatable, dimension(      :  ) :: f0_ir       ! Fraction of zeros in infrared C02 K-coefficients
  real(r8), allocatable, dimension(  :      ) :: logpint     ! Interpolated reference pressure for K-coefficients (logPa)
  real(r8), allocatable, dimension(      :  ) :: tauray_vis  ! Rayleigh scattering optical depth at each visible spectral interval
  real(r8), allocatable, dimension(      :,:) :: plnk_ir     ! Integrated Planck function lookup table at each IR spectral interval (W m-2 cm-1)
  real(r8), allocatable, dimension(      :  ) :: sol_flx     ! Solar flux within each visible spectral interval at 1AU (W m-2)

contains

  subroutine mars_nasa_rad_init(nlev)

    integer, intent(in) :: nlev

    real(r8) a, b, t
    integer i, j

    call mars_nasa_rad_final()

    nlev_rad = nlev + 1

    call mars_nasa_spectra_init()
    call mars_nasa_optics_init(nlev_rad)

    call read_kcoef()
    call interp_kcoef()

    allocate(tauray_vis(spec_vis%n     ))
    allocate(plnk_ir   (spec_ir %n,8501))
    allocate(sol_flx   (spec_vis%n     ))

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    sol_flx = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]

    do i = 1, spec_vis%n
      tauray_vis(i) = 8.7_r8 / g * (1.527_r8 * (1 + 0.013_r8 / spec_vis%wl(i)**2) / spec_vis%wl(i)**4) / (ray_p0 / 100)
    end do

    ! For each IR wavelength interval, compute the integral of B(T), the
    ! Planck function, divided by the wavelength interval, in cm-1.  The
    ! integration is in MKS units (W m^-2 cm^-1).
    do i = 1, spec_ir%n
      a = 0.5d-2 * (spec_ir%bwn(i+1) - spec_ir%bwn(i)) / (spec_ir%bwn(i) * spec_ir%bwn(i+1))
      b = 0.5d-2 * (spec_ir%bwn(i+1) + spec_ir%bwn(i)) / (spec_ir%bwn(i) * spec_ir%bwn(i+1))
      do j = 500, 9000
        t = j * 0.1_r8
        plnk_ir(i,j-499) = integrate_planck_function(a, b, t) * a / (pi * spec_ir%dwn(i))
      end do
    end do

  end subroutine mars_nasa_rad_init

  subroutine mars_nasa_rad_final()

    call mars_nasa_optics_final()

    if (allocated(tref          )) deallocate(tref          )
    if (allocated(pref          )) deallocate(pref          )
    if (allocated(qco2ref       )) deallocate(qco2ref       )
    if (allocated(qh2oref       )) deallocate(qh2oref       )
    if (allocated(gwgt          )) deallocate(gwgt          )
    if (allocated(klut_in_vis   )) deallocate(klut_in_vis   )
    if (allocated(klut_in_ir    )) deallocate(klut_in_ir    )
    if (allocated(klut_vis      )) deallocate(klut_vis      )
    if (allocated(klut_ir       )) deallocate(klut_ir       )
    if (allocated(f0_vis        )) deallocate(f0_vis        )
    if (allocated(f0_ir         )) deallocate(f0_ir         )
    if (allocated(logpint       )) deallocate(logpint       )
    if (allocated(plnk_ir       )) deallocate(plnk_ir       )
    if (allocated(sol_flx       )) deallocate(sol_flx       )
    if (allocated(tauray_vis    )) deallocate(tauray_vis    )

  end subroutine mars_nasa_rad_final

  subroutine read_kcoef()

    call fiona_open_dataset('kcoef', file_path=kcoef_file)
    call fiona_get_dim('kcoef', 'tref' , size=ntref )
    call fiona_get_dim('kcoef', 'pref' , size=npref )
    call fiona_get_dim('kcoef', 'qref' , size=nqref )
    call fiona_get_dim('kcoef', 'gauss', size=ngauss)
    allocate(tref       (ntref                             ))
    allocate(pref       (      npref                       ))
    allocate(qco2ref    (            nqref                 ))
    allocate(qh2oref    (            nqref                 ))
    allocate(gwgt       (                            ngauss))
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

  end subroutine read_kcoef

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

  pure real(r8) function integrate_planck_function(s1, s2, t) result(res)

    real(r8), intent(in) :: s1
    real(r8), intent(in) :: s2
    real(r8), intent(in) :: t   ! Temperature (K)
    
    ! C1 and C2 values from Goody and Yung (2nd edition) MKS units
    ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
    real(8), parameter :: c1 = 3.741832d-16      ! W m-2
    real(8), parameter :: c2 = 1.438786d-2       ! m K
    real(8), parameter :: x(12) = [ & ! Quadrature points
      -0.981560634246719D0,  -0.904117256370475D0, &
      -0.769902674194305D0,  -0.587317954286617D0, &
      -0.367831498998180D0,  -0.125233408511469D0, &
       0.125233408511469D0,   0.367831498998180D0, &
       0.587317954286617D0,   0.769902674194305D0, &
       0.904117256370475D0,   0.981560634246719D0  &
    ]
    real(8), parameter :: w(12) = [ & ! Quadrature weights
       0.047175336386512D0,   0.106939325995318D0, &
       0.160078328543346D0,   0.203167426723066D0, &
       0.233492536538355D0,   0.249147045813403D0, &
       0.249147045813403D0,   0.233492536538355D0, &
       0.203167426723066D0,   0.160078328543346D0, &
       0.106939325995318D0,   0.047175336386512D0  &
    ]
    real(8) wl
    integer i

    res = 0
    do i = 1, 12
      wl = s1 * x(i) + s2
      res = res + w(i) * c1 / (wl**5 * (exp(c2 / wl / t) - 1)) ! Planck function
    end do

  end function integrate_planck_function

end module mars_nasa_rad_mod
