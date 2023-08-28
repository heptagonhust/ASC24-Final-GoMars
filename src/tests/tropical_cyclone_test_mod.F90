module tropical_cyclone_test_mod

! ======================================================================
!
!  Date:  July 29, 2015
!
!  Function for setting up idealized tropical cyclone initial conditions
!
!  SUBROUTINE tropical_cyclone_sample(
!    lon,lat,p,z,zcoords,u,v,t,ptv,phis,ps,rho,q)
!
!  Given a point specified by:
!      lon    longitude (radians)
!      lat    latitude (radians)
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!
!  the functions will return:
!        p    pressure if z is specified (Pa)
!        z    geopotential height if p is specified (m)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!      ptv    virtual potential temperature (K)
!     phis    surface geopotential (m^2 s^-2)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    specific humidity (kg/kg)
!
!  Initial data are currently identical to:
!
!       Reed, K. A., and C. Jablonowski, 2011: An analytic
!       vortex initialization technique for idealized tropical
!       cyclone studies in AGCMs. Mon. Wea. Rev., 139, 689-710.
!
!  Author: Kevin A. Reed
!          Stony Brook University
!          Email: kevin.a.reed@stonybrook.edu
!
! ======================================================================
! Revision history:
!
!   2023-08: Revised by Li Dong following Yi Zhang's work.
! ======================================================================

  use flogger
  use const_mod, only: r8, pi, rad, deg, rd, g, cpd, omega, a => radius
  use namelist_mod
  use latlon_parallel_mod

  implicit none

  private

  public tropical_cyclone_test_set_diag
  public tropical_cyclone_test_set_ic

  integer , parameter :: ngauss = 20
  real(r8), parameter, dimension(ngauss), private :: gaussx = &
    [-0.0765265211334973, 0.0765265211334973, &
     -0.2277858511416451, 0.2277858511416451, &
     -0.3737060887154195, 0.3737060887154195, &
     -0.5108670019508271, 0.5108670019508271, &
     -0.6360536807265150, 0.6360536807265150, &
     -0.7463319064601508, 0.7463319064601508, &
     -0.8391169718222188, 0.8391169718222188, &
     -0.9122344282513259, 0.9122344282513259, &
     -0.9639719272779138, 0.9639719272779138, &
     -0.9931285991850949, 0.9931285991850949]
  real(r8), parameter, dimension(ngauss), private :: gaussw = &
    [0.1527533871307258 , 0.1527533871307258, &
     0.1491729864726037 , 0.1491729864726037, &
     0.1420961093183820 , 0.1420961093183820, &
     0.1316886384491766 , 0.1316886384491766, &
     0.1181945319615184 , 0.1181945319615184, &
     0.1019301198172404 , 0.1019301198172404, &
     0.0832767415767048 , 0.0832767415767048, &
     0.0626720483341091 , 0.0626720483341091, &
     0.0406014298003869 , 0.0406014298003869, &
     0.0176140071391521 , 0.0176140071391521]

!=======================================================================
!    Physical constants
!=======================================================================

  real(r8), parameter ::   &
    Lvap  = 2.5d6        , & ! Latent heat of vaporization of water
    Mvap  = 0.608d0      , & ! Ratio of molar mass of dry air/water
    p0    = 100000.0d0       ! surface pressure (Pa)

!=======================================================================
!    Test case parameters
!=======================================================================
  real(r8), parameter ::   &
    rp       = 282000.d0 , & ! Radius for calculation of PS
    dp       = 1115.d0   , & ! Delta P for calculation of PS
    zp       = 7000.d0   , & ! Height for calculation of P
    q0       = 0.021d0   , & ! q at surface from Jordan
    gamma    = 0.007d0   , & ! lapse rate
    Ts0      = 302.15d0  , & ! Surface temperature (SST)
    p00      = 101500.d0 , & ! global mean surface pressure
    cen_lat  = 10 * rad  , & ! Center latitude of initial vortex
    cen_lon  = 180 * rad , & ! Center longitufe of initial vortex
    zq1      = 3000.d0   , & ! Height 1 for q calculation
    zq2      = 8000.d0   , & ! Height 2 for q calculation
    exppr    = 1.5d0     , & ! Exponent for r dependence of p
    exppz    = 2.d0      , & ! Exponent for z dependence of p
    ztrop    = 15000.d0  , & ! Tropopause Height
    qtrop    = 1.d-11    , & ! Tropopause specific humidity
    rfpi     = 1000000.d0, & ! Radius within which to use fixed-point iteration
    constTv  = 0.608d0   , & ! Constant for Virtual Temp Conversion
    deltaz   = 2.d-13    , & ! Small number to ensure convergence in FPI
    T0       = Ts0 * (1 + constTv * q0), & ! Surface temperature
    Ttrop    = T0 - gamma * ztrop          ! Tropopause temperature
  real(r8) exponent          ! rd * gamma / g
  real(r8) ptrop             ! Tropopause pressure

contains

  subroutine tropical_cyclone_test_set_diag(blocks)

    use diag_state_mod
    use block_mod, only: block_type

    type(block_type), intent(in) :: blocks(:)

    integer iblk

    ! do iblk = 1, size(blocks)
    !   call diag_state(iblk)%init_height_levels(blocks(iblk), [100.0_r8], instance)
    ! end do

  end subroutine tropical_cyclone_test_set_diag

  real(r8) function get_dry_air_pressure(lon, lat, ptop, ztop, z) result(res)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: ptop
    real(r8), intent(in) :: ztop
    real(r8), intent(in) :: z

    integer jgw
    real(r8) z1, z2, a, b, zg
    real(r8) p, u, v, t, ptv, gzs, ps, rho, qv

    ! Set vertical height range.
    z1 = z
    z2 = ztop
    ! Transform height into Gaussian quadrature coordinate [-1,1].
    a = 0.5_r8 * (z2 - z1)
    b = 0.5_r8 * (z1 + z2)
    ! Integrate hydrostatic equation to get dry air surface pressure.
    res = 0
    do jgw = 1, ngauss
      zg = a * gaussx(jgw) + b
      call tropical_cyclone_test(lon, lat, p, zg, 1, u, v, t, ptv, gzs, ps, rho, qv)
      ! Remove water vapor from integration.
      ! Note: qv is wet mixing ratio or specific humidity here.
      res = res + gaussw(jgw) * rho * (1 - qv)
    end do
    res = a * g * res + ptop

  end function get_dry_air_pressure

  real(r8) function get_pressure(lon, lat, z) result(res)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(inout) :: z

    real(r8) u, v, t, ptv, gzs, ps, rho, qv

    call tropical_cyclone_test(lon, lat, res, z, 1, u, v, t, ptv, gzs, ps, rho, qv)

  end function get_pressure

  real(r8) function get_height(lon, lat, ptop, ztop, p) result(res)

    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: ptop
    real(r8), intent(in) :: ztop
    real(r8), intent(in) :: p    ! Dry air pressure

    real(r8), parameter :: eps = 1.0e-12_r8
    real(r8) z0, z1, z2
    real(r8) p0, p1, p2
    integer i

    z0 = 0    ; p0 = get_dry_air_pressure(lon, lat, ptop, ztop, z0)
    z1 = 10000; p1 = get_dry_air_pressure(lon, lat, ptop, ztop, z1)
    do i = 1, 1000
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
      p2 = get_dry_air_pressure(lon, lat, ptop, ztop, z2)
      if (abs(p2 - p) / p < eps .or. abs(z1 - z2) < eps .or. abs(p1 - p2) < eps) then
        exit
      end if
      z0 = z1; p0 = p1
      z1 = z2; p1 = p2
    end do
    if (i == 1001) then
      call log_error('get_height: Iteration did not converge!')
    end if
    res = z2

  end function get_height

  subroutine tropical_cyclone_test_set_ic(block)

    use block_mod
    use tracer_mod
    use formula_mod
    use vert_coord_mod

    type(block_type), intent(inout), target :: block

    integer i, j, k, jgw
    real(r8), pointer :: qv(:,:,:)
    real(r8) ztop, ps, rho, ptv

    real(r8) tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

    exponent = rd * gamma / g
    ptrop = p00 * (Ttrop / T0)**(1.d0 / exponent)

    associate (mesh   => block%mesh,             &
               lon    => block%mesh%full_lon   , &
               lat    => block%mesh%full_lat   , &
               mgs    => block%dstate(1)%mgs   , &
               mg_lev => block%dstate(1)%mg_lev, &
               ph     => block%dstate(1)%ph    , &
               ph_lev => block%dstate(1)%ph_lev, &
               z      => block%dstate(1)%gz    , &
               z_lev  => block%dstate(1)%gz_lev, &
               u      => block%dstate(1)%u     , &
               u_lon  => block%dstate(1)%u_lon , &
               v      => block%dstate(1)%v     , &
               v_lat  => block%dstate(1)%v_lat , &
               t      => block%dstate(1)%t     , &
               pt     => block%dstate(1)%pt    , &
               gzs    => block%static%gzs      )
    call tracer_get_array(block%id, idx_qv, qv, __FILE__, __LINE__)
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        ! 1. Get model top height.
        call tropical_cyclone_test(lon(i), lat(j), ptop, ztop, 0, &
          tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        ! 2. Get dry air surface pressure from temperature and moisture
        !    profiles by using Gaussian quadrature.
        mgs(i,j) = get_dry_air_pressure(lon(i), lat(j), ptop, ztop, 0.0_r8)
        ! 3. Get heights of model half levels and full levels.
        do k = mesh%half_kds, mesh%half_kde
          mg_lev(i,j,k) = vert_coord_calc_mg_lev(k, mgs(i,j))
          z_lev(i,j,k) = get_height(lon(i), lat(j), ptop, ztop, mg_lev(i,j,k))
        end do
        do k = mesh%full_kds, mesh%full_kde
          z(i,j,k) = 0.5_r8 * (z_lev(i,j,k) + z_lev(i,j,k+1))
        end do
        ! 4. Get variables on model full levels.
        do k = mesh%full_kds, mesh%full_kde
          call tropical_cyclone_test(lon(i), lat(j), ph(i,j,k), z(i,j,k), 1, &
            u(i,j,k), v(i,j,k), t(i,j,k), ptv, gzs(i,j), ps, rho, qv(i,j,k))
          ! Convert to dry mixing ratio.
          qv(i,j,k) = qv(i,j,k) / (1 - qv(i,j,k))
          pt(i,j,k) = modified_potential_temperature(t(i,j,k), ph(i,j,k), qv(i,j,k))
        end do
      end do
    end do
    call fill_halo(block%halo,         u, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block%halo,         v, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
    call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          u_lon(i,j,k) = 0.5_r8 * (u(i,j,k) + u(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(block%halo, u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat(i,j,k) = 0.5_r8 * (v(i,j,k) + v(i,j+1,k))
        end do
      end do
    end do
    call fill_halo(block%halo, v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
    if (pole_damp_mgs) then
      call fill_halo(block%filter_halo, mgs, full_lon=.true., full_lat=.true.)
    else
      call fill_halo(block%halo, mgs, full_lon=.true., full_lat=.true.)
    end if
    call tracer_calc_qm(block)
    end associate

  end subroutine tropical_cyclone_test_set_ic

!=======================================================================
!    Evaluate the tropical cyclone initial conditions
!=======================================================================
  subroutine tropical_cyclone_test(lon, lat, p, z, zcoords, u, v, t, &
                                   ptv, phis, ps, rho, q)

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------

    real(r8), intent(in) ::    &
              lon,             &     ! Longitude (radians)
              lat                    ! Latitude (radians)

    real(r8), intent(inout) :: &
              p,               &     ! Pressure (Pa)
              z                      ! Height (m)

    integer, intent(in) :: zcoords   ! 1 if z coordinates are specified
                                     ! 0 if p coordinates are specified

    real(r8), intent(out) ::   &
              u,               &     ! Zonal wind (m s^-1)
              v,               &     ! Meridional wind (m s^-1)
              t,               &     ! Temperature (K)
              ptv,             &     ! Virtual potential temperature (K)
              phis,            &     ! Surface Geopotential (m^2 s^-2)
              ps,              &     ! Surface Pressure (Pa)
              rho,             &     ! Density (kg m^-3)
              q                      ! Specific Humidity (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    real(r8)  d1, d2, d, vfac, ufac, gr, f, zn

    integer   n

    !------------------------------------------------
    !   Define Great circle distance (gr) and
    !   Coriolis parameter (f)
    !------------------------------------------------
    ! Coriolis parameter
    f  = 2 * omega * sin(cen_lat)
    ! Great circle radius
    gr = a * acos(sin(cen_lat) * sin(lat) + &
         (cos(cen_lat) * cos(lat) * cos(lon - cen_lon)))

    !------------------------------------------------
    !   Initialize PS (surface pressure)
    !------------------------------------------------
    ps = p00 - dp * exp(-(gr / rp)**exppr)

    !------------------------------------------------
    !   Initialize altitude (z) if pressure provided
    !   or pressure if altitude (z) is provided
    !------------------------------------------------
    if (zcoords == 1) then
      ! Eqn (8):
      if (z > ztrop) then
        p = ptrop * exp(-(g * (z - ztrop))/(rd * Ttrop))
      else
        p = (p00 - dp * exp(-(gr / rp)**exppr) &
            * exp(-(z / zp)**exppz)) &
            * ((T0 - gamma * z) / T0)**(1 / exponent)
      end if
    else
      ! Eqn (24):
      if (ps >= p .and. p >= ptrop) then
        z = (T0 / gamma) * (1 - (p / ps)**exponent)
      else if (ptrop >= p) then
        z = ztrop + rd * Ttrop / g * log(ptrop / p)
      else
        call log_error('p < ps!', __FILE__, __LINE__)
      end if
      ! If inside a certain distance of the center of the storm
      ! perform a Fixed-point iteration to calculate the height
      ! more accurately
      if (gr < rfpi) then
        do n = 1, 21
          zn = z - fpif(p, gr, z) / fpidfdz(gr, z)
          if (n == 21) then
            write(*, *) 'FPI did not converge after 20 interations in q & T!!!'
          else if (abs(zn - z) / abs(zn) > deltaz) then
            z = zn
          else
            exit
          end if
        end do
        z = zn

        ! n = 1
        ! 20 continue
        ! n = n + 1
        ! zn = z - fpif(p, gr, z) / fpidfdz(gr, z)
        ! if (n > 20) then
        !   write(*, *) 'FPI did not converge after 20 interations in q & T!!!'
        ! else if (abs(zn - z) / abs(zn) > deltaz) then
        !   z = zn
        !   goto 20
        ! end if
        ! z = zn
      end if
    end if

    !------------------------------------------------
    !   Initialize U and V (wind components)
    !------------------------------------------------
    d1 = sin(cen_lat) * cos(lat) - &
         cos(cen_lat) * sin(lat) * cos(lon - cen_lon)
    d2 = cos(cen_lat) * sin(lon - cen_lon)
    d  = max(1.0e-25_r8, sqrt(d1**2 + d2**2))
    ufac = d1 / d
    vfac = d2 / d

    if (z > ztrop) then
      u = 0
      v = 0
    else
      v = vfac * (-f * gr / 2.d0 + sqrt((f * gr / 2.d0)**2 &
          - exppr * (gr / rp)**exppr * rd * (T0 - gamma * z) &
          / (exppz * z * rd * (T0 - gamma * z) / (g * zp**exppz) &
          + (1 - p00 / dp * exp((gr / rp)**exppr) * exp((z / zp)**exppz)))))
      u = ufac * (-f * gr / 2.d0 + sqrt((f * gr / 2.d0)**2 &
          - exppr * (gr / rp)**exppr * rd * (T0 - gamma * z) &
          / (exppz * z * rd * (T0 - gamma * z) / (g * zp**exppz) &
          + (1 - p00 / dp * exp((gr / rp)**exppr) * exp((z / zp)**exppz)))))
    end if

    !------------------------------------------------
    !   Initialize water vapor mixing ratio (q)
    !------------------------------------------------
    if (z > ztrop) then
      q = qtrop
    else
      q = q0 * exp(-z / zq1) * exp(-(z / zq2)**exppz)
    end if

    !------------------------------------------------
    !   Initialize temperature (T)
    !------------------------------------------------
    if (z > ztrop) then
      t = Ttrop
    else
      t = (T0 - gamma * z) / (1 + constTv * q) &
          / (1 + exppz * rd * (T0 - gamma * z) * z &
          / (g * zp**exppz * (1 - p00 / dp * exp((gr / rp)**exppr) &
          * exp((z / zp)**exppz))))
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature (ptv)
    !-----------------------------------------------------
    ptv = t * (1 + constTv * q) * (p0 / p)**(rd / cpd)

    !-----------------------------------------------------
    !   Initialize surface geopotential (PHIS)
    !-----------------------------------------------------
    phis = 0

    !-----------------------------------------------------
    !   Initialize density (rho)
    !-----------------------------------------------------
    rho = p / (rd * t * (1 + constTv * q))

  end subroutine tropical_cyclone_test

!-----------------------------------------------------------------------
!    First function for fixed point iterations
!-----------------------------------------------------------------------
  real(r8) function fpif(p, gr, z) result(res)

    real(r8), intent(in) :: p, gr, z

    if (z >= 0 .and. z <= ztrop) then
      res = p - (p00 - dp * exp(-(gr / rp)**exppr) &
          * exp(-(z / zp)**exppz)) &
          * ((T0 - gamma * z) / T0)**(g / (rd * gamma))
    else
      res = p - (ptrop * exp(-g * (z - ztrop) / (rd * Ttrop)))
    end if

  end function fpif

!-----------------------------------------------------------------------
!    Second function for fixed point iterations
!-----------------------------------------------------------------------
  real(r8) function fpidfdz(gr, z) result(res)

    real(r8), intent(in) :: gr, z

    res = -exppz * dp * z / (zp**2) * exp(-(gr / rp)**exppr) &
        * exp(-(z / zp)**exppz) * ((T0 - gamma * z) / T0)**(g / (rd * gamma)) &
        + g / (rd * T0) * (p00 - dp * exp(-(gr / rp)**exppr) * exp(-(z / zp)**exppz)) &
        * ((T0 - gamma * z) / T0)**(g / (rd * gamma) - 1)

  end function fpidfdz

end module tropical_cyclone_test_mod
