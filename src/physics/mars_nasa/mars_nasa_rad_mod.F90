module mars_nasa_rad_mod

  use fiona
  use const_mod
  use mars_nasa_namelist_mod

  implicit none

  integer               nlev_rad
  integer, parameter :: nspec_vis = 7  ! Number of visible spectral intervals
  integer, parameter :: nspec_ir  = 5  ! Number of IR spectral intervals
  integer               ntref          ! Number of reference temperature levels for K-coefficients
  integer               npref          ! Number of reference pressure levels for K-coefficients
  integer               nqref          ! Number of different water mixing ratios for K-coefficients that are now CO2 + H2O
  integer               ngauss         ! Number of Gaussian quadrature points for K-coefficients
  integer               npint          ! Number of Lagrange interpolated reference pressures for the CO2 K-coefficients
  integer               nbin_dust      ! Number of dust particle size bins
  integer, parameter :: taumax    = 35 ! Maximum optical depth

  ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
  ! spectral bins.  Go from smaller to larger wavenumbers, the same as
  ! in the IR.
  real(r8), parameter :: bwn_vis(nspec_vis+1) = [ &
     2222.22_r8, & ! ->   4.50 microns
     3087.37_r8, & ! ->   3.24 microns
     4030.63_r8, & ! ->   2.48 microns
     5370.57_r8, & ! ->   1.86 microns
     7651.11_r8, & ! ->   1.31 microns
    12500.00_r8, & ! ->   0.80 microns
    25000.00_r8, & ! ->   0.40 microns
    41666.67_r8  & ! ->   0.24 microns
  ]
  ! Bin wavenumber of the edges of the IR spectral bins
  ! units are inverse centimeters.  Dimension needs to be changed
  ! if the number of IR bins changes.
  real(r8), parameter :: bwn_ir(nspec_ir+1) = [ &
      10.000_r8, & ! -> 1000.00 microns
     166.667_r8, & ! ->   60.00 microns
     416.667_r8, & ! ->   24.00 microns
     833.333_r8, & ! ->   12.00 microns
    1250.000_r8, & ! ->    8.00 microns
    2222.222_r8  & ! ->    4.50 microns
  ]
  real(r8), parameter :: ray_p0   = 9.423e6_r8        ! Reference pressure for Rayleigh scattering (Pa)

  ! Dust optical properties for log-normal size distribution with reff = 1.5 micron, veff = 0.5
  real(r8), parameter :: dust_qext0_vis(nspec_vis) = [ &
    1.834_r8, &
    2.296_r8, &
    2.672_r8, &
    2.829_r8, &
    2.698_r8, &
    2.452_r8, &
    2.261_r8  &
  ]
  real(r8), parameter :: dust_qscat0_vis(nspec_vis) = [ &
    1.695_r8, &
    2.031_r8, &
    2.583_r8, &
    2.744_r8, &
    2.626_r8, &
    2.225_r8, &
    1.525_r8  &
  ]
  real(r8), parameter :: dust_g0_vis(nspec_vis) = [ &
    0.551_r8, &
    0.640_r8, &
    0.661_r8, &
    0.678_r8, &
    0.690_r8, &
    0.743_r8, &
    0.868_r8  &
  ]
  real(r8), parameter :: dust_qext0_ir(nspec_ir) = [ &
    0.008_r8, &
    0.262_r8, &
    0.491_r8, &
    1.017_r8, &
    0.444_r8  &
  ]
  real(r8), parameter :: dust_qscat0_ir(nspec_ir) = [ &
    0.001_r8, &
    0.037_r8, &
    0.122_r8, &
    0.351_r8, &
    0.336_r8  &
  ]
  real(r8), parameter :: dust_g0_ir(nspec_ir) = [ &
    0.004_r8, &
    0.030_r8, &
    0.095_r8, &
    0.214_r8, &
    0.316_r8  &
  ]
  real(r8), allocatable, dimension(:,:) :: dust_qext_vis
  real(r8), allocatable, dimension(:,:) :: dust_qscat_vis
  real(r8), allocatable, dimension(:,:) :: dust_g_vis
  real(r8), allocatable, dimension(:,:) :: dust_qext_ir
  real(r8), allocatable, dimension(:,:) :: dust_qscat_ir
  real(r8), allocatable, dimension(:,:) :: dust_g_ir

  real(r8), allocatable, dimension(:        ) :: tref         ! Reference temperature for K-coefficients (K)
  real(r8), allocatable, dimension(  :      ) :: pref         ! Reference pressure for K-coefficients (Pa)
  real(r8), allocatable, dimension(    :    ) :: qh2oref      ! Reference water vapor volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(    :    ) :: qco2ref      ! Reference CO2 volume mixing ratio for K-coefficients (1)
  real(r8), allocatable, dimension(        :) :: gwgt         ! Gaussian point weights
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_vis  ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_in_ir   ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_vis     ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: klut_ir      ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(      :  ) :: f0_vis       ! Fraction of zeros in visible CO2 K-coefficients
  real(r8), allocatable, dimension(      :  ) :: f0_ir        ! Fraction of zeros in infrared C02 K-coefficients
  real(r8), allocatable, dimension(  :      ) :: logpint      ! Interpolated reference pressure for K-coefficients (logPa)

  real(r8), allocatable, dimension(      :  ) :: wn_vis       ! Wavenumbers at spectral interval centers for visible (cm-1)
  real(r8), allocatable, dimension(      :  ) :: wn_ir        ! Wavenumbers at spectral interval centers for IR (cm-1)
  real(r8), allocatable, dimension(      :  ) :: dwn_vis      ! Delta wavenumber of each visible spectral interval (cm-1)
  real(r8), allocatable, dimension(      :  ) :: dwn_ir       ! Delta wavenumber of each IR spectral interval (cm-1)
  real(r8), allocatable, dimension(      :  ) :: wl_vis       ! Wavelengths at spectral interval centers for visible (micron)
  real(r8), allocatable, dimension(      :  ) :: wl_ir        ! Wavelengths at spectral interval centers for IR (micron)
  real(r8), allocatable, dimension(      :  ) :: qext_vis     ! Extinction efficiency at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: qext_ir      ! Extinction efficiency at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: qscat_vis    ! Scattering efficiency at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: qscat_ir     ! Scattering efficiency at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: wscat_vis    ! Single scattering albedo at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: wscat_ir     ! Single scattering albedo at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: gscat_vis    ! Asymmetry parameter at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: gscat_ir     ! Asymmetry parameter at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: tauray_vis   ! Rayleigh scattering optical depth at each visible spectral interval
  real(r8), allocatable, dimension(      :,:) :: plnk_ir      ! Integrated Planck function lookup table at each IR spectral interval (W m-2 cm-1)
  real(r8), allocatable, dimension(      :  ) :: sol_flx      ! Solar flux within each visible spectral interval at 1AU (W m-2)

contains

  subroutine mars_nasa_rad_init(nlev)

    integer, intent(in) :: nlev

    real(r8) a, b, t
    integer i, j

    call mars_nasa_rad_final()

    nlev_rad = nlev + 1

    allocate(wn_vis     (                  nspec_vis       ))
    allocate(wn_ir      (                  nspec_ir        ))
    allocate(dwn_vis    (                  nspec_vis       ))
    allocate(dwn_ir     (                  nspec_ir        ))
    allocate(wl_vis     (                  nspec_vis       ))
    allocate(wl_ir      (                  nspec_ir        ))
    allocate(qext_vis   (                  nspec_vis       ))
    allocate(qext_ir    (                  nspec_ir        ))
    allocate(qscat_vis  (                  nspec_vis       ))
    allocate(qscat_ir   (                  nspec_ir        ))
    allocate(wscat_vis  (                  nspec_vis       ))
    allocate(wscat_ir   (                  nspec_ir        ))
    allocate(gscat_vis  (                  nspec_vis       ))
    allocate(gscat_ir   (                  nspec_ir        ))
    allocate(tauray_vis (                  nspec_vis       ))
    allocate(plnk_ir    (                  nspec_ir ,8501  ))
    allocate(sol_flx    (                  nspec_vis       ))

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    sol_flx = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]

    ! Visible spectral intervals
    do i = 1, nspec_vis
      wn_vis (i) = 0.5_r8 * (bwn_vis(i) + bwn_vis(i+1))
      dwn_vis(i) = bwn_vis(i+1) - bwn_vis(i)
      wl_vis (i) = 1.0e4_r8 / wn_vis(i)
      tauray_vis(i) = 8.7_r8 / g * (1.527_r8 * (1 + 0.013_r8 / wl_vis(i)**2) / wl_vis(i)**4) / (ray_p0 / 100)
    end do

    ! IR spectral intervals
    do i = 1, nspec_ir
      wn_ir (i) = 0.5_r8 * (bwn_ir(i) + bwn_ir(i+1))
      dwn_ir(i) = bwn_ir(i+1) - bwn_ir(i)
      wl_ir (i) = 1.0e4_r8 / wn_ir(i)
    end do

    ! For each IR wavelength interval, compute the integral of B(T), the
    ! Planck function, divided by the wavelength interval, in cm-1.  The
    ! integration is in MKS units (W m^-2 cm^-1).
    do i = 1, nspec_ir
      a = 0.5d-2 * (bwn_ir(i+1) - bwn_ir(i)) / (bwn_ir(i) * bwn_ir(i+1))
      b = 0.5d-2 * (bwn_ir(i+1) + bwn_ir(i)) / (bwn_ir(i) * bwn_ir(i+1))
      do j = 500, 9000
        t = j * 0.1_r8
        plnk_ir(i,j-499) = integrate_planck_function(a, b, t) * a / (pi * dwn_ir(i))
      end do
    end do

    call read_kcoef()
    call interp_kcoef()
    call read_dust_optics()

    do i = 1, nspec_vis
      qext_vis(i) = dust_qext0_vis(i)
      qscat_vis(i) = dust_qscat0_vis(i)
      if (qscat_vis(i) >= qext_vis(i)) then
        qscat_vis(i) = 0.99999_r8 * qext_vis(i)
      end if
      wscat_vis(i) = qscat_vis(i) / qext_vis(i)
      gscat_vis(i) = dust_g0_vis(i)
    end do

    do i = 1, nspec_ir
      qext_ir(i) = dust_qext0_ir(i)
      qscat_ir(i) = dust_qscat0_ir(i)
      if (qscat_ir(i) >= qext_ir(i)) then
        qscat_ir(i) = 0.99999_r8 * qext_ir(i)
      end if
      wscat_ir(i) = qscat_ir(i) / qext_ir(i)
      gscat_ir(i) = dust_g0_ir(i)
    end do

  end subroutine mars_nasa_rad_init

  subroutine mars_nasa_rad_final()

    if (allocated(dust_qext_vis )) deallocate(dust_qext_vis )
    if (allocated(dust_qscat_vis)) deallocate(dust_qscat_vis)
    if (allocated(dust_g_vis    )) deallocate(dust_g_vis    )
    if (allocated(dust_qext_ir  )) deallocate(dust_qext_ir  )
    if (allocated(dust_qscat_ir )) deallocate(dust_qscat_ir )
    if (allocated(dust_g_ir     )) deallocate(dust_g_ir     )
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
    if (allocated(wn_vis        )) deallocate(wn_vis        )
    if (allocated(wn_ir         )) deallocate(wn_ir         )
    if (allocated(dwn_vis       )) deallocate(dwn_vis       )
    if (allocated(dwn_ir        )) deallocate(dwn_ir        )
    if (allocated(wl_vis        )) deallocate(wl_vis        )
    if (allocated(wl_ir         )) deallocate(wl_ir         )
    if (allocated(qext_vis      )) deallocate(qext_vis      )
    if (allocated(qext_ir       )) deallocate(qext_ir       )
    if (allocated(qscat_vis     )) deallocate(qscat_vis     )
    if (allocated(qscat_ir      )) deallocate(qscat_ir      )
    if (allocated(wscat_vis     )) deallocate(wscat_vis     )
    if (allocated(wscat_ir      )) deallocate(wscat_ir      )
    if (allocated(gscat_vis     )) deallocate(gscat_vis     )
    if (allocated(gscat_ir      )) deallocate(gscat_ir      )
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
    allocate(klut_in_vis(ntref,npref,nqref,nspec_vis,ngauss))
    allocate(klut_in_ir (ntref,npref,nqref,nspec_ir ,ngauss))
    allocate(f0_vis     (                  nspec_vis       ))
    allocate(f0_ir      (                  nspec_ir        ))
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

  subroutine read_dust_optics()

    call fiona_open_dataset('dust_optics', file_path=dust_optics_file)
    call fiona_get_dim('dust_optics', 'bin', size=nbin_dust)
    allocate(dust_qext_vis (nbin_dust,nspec_vis))
    allocate(dust_qscat_vis(nbin_dust,nspec_vis))
    allocate(dust_g_vis    (nbin_dust,nspec_vis))
    allocate(dust_qext_ir  (nbin_dust,nspec_ir ))
    allocate(dust_qscat_ir (nbin_dust,nspec_ir ))
    allocate(dust_g_ir     (nbin_dust,nspec_ir ))
    call fiona_start_input('dust_optics')
    call fiona_input('dust_optics', 'dust_qext_vis' , dust_qext_vis )
    call fiona_input('dust_optics', 'dust_qscat_vis', dust_qscat_vis)
    call fiona_input('dust_optics', 'dust_g_vis'    , dust_g_vis    )
    call fiona_input('dust_optics', 'dust_qext_ir'  , dust_qext_ir  )
    call fiona_input('dust_optics', 'dust_qscat_ir' , dust_qscat_ir )
    call fiona_input('dust_optics', 'dust_g_ir'     , dust_g_ir     )
    call fiona_end_input('dust_optics')

  end subroutine read_dust_optics

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
      do iw = 1, nspec_vis
        if (klut_in_vis(it,ip,iq,iw,ig) > 1.0d-200) then
          klut_in_vis(it,ip,iq,iw,ig) = log10(klut_in_vis(it,ip,iq,iw,ig))
        else
          klut_in_vis(it,ip,iq,iw,ig) = -200
        end if
      end do
      do iw = 1, nspec_ir
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

    allocate(klut_vis(ntref,npint,nqref,nspec_vis,ngauss))
    allocate(klut_ir (ntref,npint,nqref,nspec_ir ,ngauss))

    do ig = 1, ngauss
    do iw = 1, nspec_vis
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
    do iw = 1, nspec_ir
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
