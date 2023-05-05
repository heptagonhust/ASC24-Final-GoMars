module mars_rad_nasa_mod

  use fiona
  use const_mod
  use namelist_mod
  use process_mod

  implicit none

  integer, parameter :: nspec_vis = 7                 ! Number of visible spectral intervals
  integer, parameter :: nspec_ir  = 5                 ! Number of IR spectral intervals
  integer               ntref                         ! Number of reference temperature levels for K-coefficients
  integer               npref                         ! Number of reference pressure levels for K-coefficients
  integer               nqref                         ! Number of different water mixing ratios for K-coefficients that are now CO2 + H2O
  integer               ngauss                        ! Number of Gaussian quadrature points for K-coefficients
  integer, parameter :: npint     = 51                ! Number of Lagrange interpolated reference pressures for the CO2 K-coefficients
  integer, parameter :: taumax    = 35                ! Maximum optical depth

  ! Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
  ! spectral bins.  Go from smaller to larger wavenumbers, the same as
  ! in the IR.
  real(r8), parameter :: bwnv(nspec_vis+1) = [ &
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
  real(r8), parameter :: bwni(nspec_ir+1) = [ &
      10.000_r8, & ! -> 1000.00 microns
     166.667_r8, & ! ->   60.00 microns
     416.667_r8, & ! ->   24.00 microns
     833.333_r8, & ! ->   12.00 microns
    1250.000_r8, & ! ->    8.00 microns
    2222.222_r8  & ! ->    4.50 microns
  ]
  real(r8), parameter :: ray_p0   = 9.423e6_r8        ! Reference pressure for Rayleigh scattering (Pa)

  ! Dust optical properties for log-normal size distribution with reff = 1.5 micron, veff = 0.5
  real(r8), parameter :: vis_dust_qext(nspec_vis) = [ &
    1.834_r8, &
    2.296_r8, &
    2.672_r8, &
    2.829_r8, &
    2.698_r8, &
    2.452_r8, &
    2.261_r8  &
  ]
  real(r8), parameter :: vis_dust_qscat(nspec_vis) = [ &
    1.695_r8, &
    2.031_r8, &
    2.583_r8, &
    2.744_r8, &
    2.626_r8, &
    2.225_r8, &
    1.525_r8  &
  ]
  real(r8), parameter :: vis_dust_g(nspec_vis) = [ &
    0.551_r8, &
    0.640_r8, &
    0.661_r8, &
    0.678_r8, &
    0.690_r8, &
    0.743_r8, &
    0.868_r8  &
  ]
  real(r8), parameter :: ir_dust_qext(nspec_ir) = [ &
    0.008_r8, &
    0.262_r8, &
    0.491_r8, &
    1.017_r8, &
    0.444_r8  &
  ]
  real(r8), parameter :: ir_dust_qscat(nspec_ir) = [ &
    0.001_r8, &
    0.037_r8, &
    0.122_r8, &
    0.351_r8, &
    0.336_r8  &
  ]
  real(r8), parameter :: ir_dust_g(nspec_ir) = [ &
    0.004_r8, &
    0.030_r8, &
    0.095_r8, &
    0.214_r8, &
    0.316_r8  &
  ]

  real(r8), allocatable, dimension(      :  ) :: sol_flx      ! Solar flux within each visible spectral interval at 1AU (W m-2)
  real(r8), allocatable, dimension(      :  ) :: vis_wn       ! Wavenumbers at spectral interval centers for visible (cm-1)
  real(r8), allocatable, dimension(      :  ) :: vis_dwn      ! Delta wavenumber of each visible spectral interval (cm-1)
  real(r8), allocatable, dimension(      :  ) :: vis_wl       ! Wavelengths at spectral interval centers for visible (micron)
  real(r8), allocatable, dimension(      :  ) :: vis_ray_tau  ! Rayleigh scattering optical depth at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: vis_qext     ! Extinction efficiency at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: vis_qscat    ! Scattering efficiency at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: vis_scat_w   ! Single scattering albedo at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: vis_scat_g   ! Asymmetry parameter at each visible spectral interval
  real(r8), allocatable, dimension(      :  ) :: ir_wn        ! Wavenumbers at spectral interval centers for IR (cm-1)
  real(r8), allocatable, dimension(      :  ) :: ir_dwn       ! Delta wavenumber of each IR spectral interval (cm-1)
  real(r8), allocatable, dimension(      :  ) :: ir_wl        ! Wavelengths at spectral interval centers for IR (micron)
  real(r8), allocatable, dimension(      :,:) :: ir_plnk      ! Planck function at each IR spectral interval (W m-2 cm-1)
  real(r8), allocatable, dimension(      :  ) :: ir_qext      ! Extinction efficiency at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: ir_qscat     ! Scattering efficiency at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: ir_scat_w    ! Single scattering albedo at each IR spectral interval
  real(r8), allocatable, dimension(      :  ) :: ir_scat_g    ! Asymmetry parameter at each IR spectral interval
  ! K-coefficients
  real(r8), allocatable, dimension(:        ) :: tref         ! Reference temperature for K-coefficients (K)
  real(r8), allocatable, dimension(  :      ) :: pref         ! Reference pressure for K-coefficients (Pa)
  real(r8), allocatable, dimension(    :    ) :: qref         ! Reference water vapor volume mixing ratio for K-coefficeints (1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: vis_k_co2    ! CO2 K-coefficients for each visible spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(:,:,:,:,:) :: ir_k_co2     ! CO2 K-coefficients for each IR spectral interval (cm2 mole-1)
  real(r8), allocatable, dimension(        :) :: gauss_wgt    ! Gaussian point weights

contains

  subroutine mars_rad_nasa_init()

    real(r8) a, b, t
    integer i, j

    call mars_rad_nasa_final()

    allocate(sol_flx    (                  nspec_vis       ))
    allocate(vis_wn     (                  nspec_vis       ))
    allocate(vis_dwn    (                  nspec_vis       ))
    allocate(vis_wl     (                  nspec_vis       ))
    allocate(vis_ray_tau(                  nspec_vis       ))
    allocate(vis_qext   (                  nspec_vis       ))
    allocate(vis_qscat  (                  nspec_vis       ))
    allocate(vis_scat_w (                  nspec_vis       ))
    allocate(vis_scat_g (                  nspec_vis       ))
    allocate(ir_wn      (                  nspec_ir        ))
    allocate(ir_dwn     (                  nspec_ir        ))
    allocate(ir_wl      (                  nspec_ir        ))
    allocate(ir_plnk    (                  nspec_ir ,8501  ))
    allocate(ir_qext    (                  nspec_ir        ))
    allocate(ir_qscat   (                  nspec_ir        ))
    allocate(ir_scat_w  (                  nspec_ir        ))
    allocate(ir_scat_g  (                  nspec_ir        ))

    call mars_rad_nasa_read_kcoef()

    ! Sum equals 1356 W m-2 (values from Wehrli, 1985)
    sol_flx = [12.7_r8, 24.2_r8, 54.6_r8, 145.9_r8, 354.9_r8, 657.5_r8, 106.3_r8]

    ! Visible spectral intervals
    do i = 1, nspec_vis
      vis_wn     (i) = 0.5_r8 * (bwnv(i) + bwnv(i+1))
      vis_dwn    (i) = bwnv(i+1) - bwnv(i)
      vis_wl     (i) = 1.0e4_r8 / vis_wn(i)
      vis_ray_tau(i) = 8.7_r8 / g * (1.527_r8 * (1 + 0.013_r8 / vis_wl(i)**2) / vis_wl(i)**4) * 100 / ray_p0
    end do

    ! IR spectral intervals
    do i = 1, nspec_ir
      ir_wn (i) = 0.5_r8 * (bwni(i) + bwni(i+1))
      ir_dwn(i) = bwni(i+1) - bwni(i)
      ir_wl (i) = 1.0e4_r8 / ir_wn(i)
    end do

    ! For each IR wavelength interval, compute the integral of B(T), the
    ! Planck function, divided by the wavelength interval, in cm-1.  The
    ! integration is in MKS units (W m^-2 cm^-1).
    do i = 1, nspec_ir
      a = 0.5d-2 * (bwni(i+1) - bwni(i)) / (bwni(i) * bwni(i+1))
      b = 0.5d-2 * (bwni(i+1) + bwni(i)) / (bwni(i) * bwni(i+1))
      do j = 500, 9000
        t = j * 0.1_r8
        ir_plnk(i,j+499) = integrate_planck_function(a, b, t) * a / (pi * ir_dwn(i))
      end do
    end do

    do i = 1, nspec_vis
      vis_qext(i) = vis_dust_qext(i)
      vis_qscat(i) = vis_dust_qscat(i)
      if (vis_qscat(i) >= vis_qext(i)) then
        vis_qscat(i) = 0.99999_r8 * vis_qext(i)
      end if
      vis_scat_w(i) = vis_qscat(i) / vis_qext(i)
      vis_scat_g(i) = vis_dust_g(i)
    end do

    do i = 1, nspec_ir
      ir_qext(i) = ir_dust_qext(i)
      ir_qscat(i) = ir_dust_qscat(i)
      if (ir_qscat(i) >= ir_qext(i)) then
        ir_qscat(i) = 0.99999_r8 * ir_qext(i)
      end if
      ir_scat_w(i) = ir_qscat(i) / ir_qext(i)
      ir_scat_g(i) = ir_dust_g(i)
    end do

  end subroutine mars_rad_nasa_init

  subroutine mars_rad_nasa_final()

    if (allocated(sol_flx    )) deallocate(sol_flx    )
    if (allocated(vis_wn     )) deallocate(vis_wn     )
    if (allocated(vis_dwn    )) deallocate(vis_dwn    )
    if (allocated(vis_wl     )) deallocate(vis_wl     )
    if (allocated(vis_ray_tau)) deallocate(vis_ray_tau)
    if (allocated(vis_qext   )) deallocate(vis_qext   )
    if (allocated(vis_qscat  )) deallocate(vis_qscat  )
    if (allocated(vis_scat_w )) deallocate(vis_scat_w )
    if (allocated(vis_scat_g )) deallocate(vis_scat_g )
    if (allocated(ir_wn      )) deallocate(ir_wn      )
    if (allocated(ir_dwn     )) deallocate(ir_dwn     )
    if (allocated(ir_wl      )) deallocate(ir_wl      )
    if (allocated(ir_plnk    )) deallocate(ir_plnk    )
    if (allocated(ir_qext    )) deallocate(ir_qext    )
    if (allocated(ir_qscat   )) deallocate(ir_qscat   )
    if (allocated(ir_scat_w  )) deallocate(ir_scat_w  )
    if (allocated(ir_scat_g  )) deallocate(ir_scat_g  )
    if (allocated(tref       )) deallocate(tref       )
    if (allocated(pref       )) deallocate(pref       )
    if (allocated(qref       )) deallocate(qref       )
    if (allocated(vis_k_co2  )) deallocate(vis_k_co2  )
    if (allocated(ir_k_co2   )) deallocate(ir_k_co2   )
    if (allocated(gauss_wgt  )) deallocate(gauss_wgt  )

  end subroutine mars_rad_nasa_final

  subroutine mars_rad_nasa_read_kcoef()

    real(8), allocatable, dimension(:,:,:,:,:) :: kv, ki

    call fiona_open_dataset('kcoef', trim(gmcore_data_dir) // '/nasa_kcoef.nc', mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_get_dim('kcoef', 't', size=ntref)
    call fiona_get_dim('kcoef', 'p', size=npref)
    call fiona_get_dim('kcoef', 'qv', size=nqref)
    call fiona_get_dim('kcoef', 'gp', size=ngauss)
    allocate(tref     (ntref                             ))
    allocate(pref     (      npref                       ))
    allocate(qref     (            nqref                 ))
    allocate(kv       (ntref,npref,nqref,nspec_vis,ngauss))
    allocate(ki       (ntref,npref,nqref,nspec_ir ,ngauss))
    allocate(vis_k_co2(ntref,npint,nqref,nspec_vis,ngauss))
    allocate(ir_k_co2 (ntref,npint,nqref,nspec_ir ,ngauss))
    allocate(gauss_wgt(                            ngauss))
    call fiona_start_input('kcoef')
    call fiona_input('kcoef', 't' , tref     )
    call fiona_input('kcoef', 'p' , pref     )
    call fiona_input('kcoef', 'qv', qref     )
    call fiona_input('kcoef', 'kv', kv       )
    call fiona_input('kcoef', 'ki', ki       )
    call fiona_input('kcoef', 'gw', gauss_wgt)
    call fiona_close_dataset('kcoef')

    deallocate(kv, ki)

    ! Interpolate K-coefficients to the refined pressure values.

  end subroutine mars_rad_nasa_read_kcoef

  pure real(r8) function integrate_planck_function(s1, s2, t) result(res)

    real(r8), intent(in) :: s1
    real(r8), intent(in) :: s2
    real(r8), intent(in) :: t
    
    ! C1 and C2 values from Goody and Yung (2nd edition) MKS units
    ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
    real(8), parameter :: c1 = 3.741832d-16      ! W m-2
    real(8), parameter :: c2 = 1.438786d-2       ! m K
    real(8), parameter :: x(12) = [ &
      -0.981560634246719D0,  -0.904117256370475D0, &
      -0.769902674194305D0,  -0.587317954286617D0, &
      -0.367831498998180D0,  -0.125233408511469D0, &
       0.125233408511469D0,   0.367831498998180D0, &
       0.587317954286617D0,   0.769902674194305D0, &
       0.904117256370475D0,   0.981560634246719D0  &
    ]
    real(8), parameter :: w(12) = [ &
       0.047175336386512D0,   0.106939325995318D0, &
       0.160078328543346D0,   0.203167426723066D0, &
       0.233492536538355D0,   0.249147045813403D0, &
       0.249147045813403D0,   0.233492536538355D0, &
       0.203167426723066D0,   0.160078328543346D0, &
       0.106939325995318D0,   0.047175336386512D0  &
    ]
    real(8) y
    integer i

    res = 0
    do i = 1, 12
      y = s1 * x(i) + s2
      res = res + w(i) * c1 / (y**5 * (exp(c2 / y / t) - 1))
    end do

  end function integrate_planck_function

end module mars_rad_nasa_mod