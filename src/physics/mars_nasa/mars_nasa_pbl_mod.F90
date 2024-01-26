module mars_nasa_pbl_mod

  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_objects_mod
  use mars_nasa_sfc_mod

  implicit none

  private

  public mars_nasa_pbl_init
  public mars_nasa_pbl_final

contains

  subroutine mars_nasa_pbl_init()

    integer iblk

    do iblk = 1, size(objects)
      associate (state => objects(iblk)%state)
      state%saved_shear2 = 0
      end associate
    end do

  end subroutine mars_nasa_pbl_init

  subroutine mars_nasa_pbl_final()

  end subroutine mars_nasa_pbl_final

  subroutine mars_nasa_pbl_run(state, tend, dt)

    type(mars_nasa_state_type), intent(inout) :: state
    type(mars_nasa_tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    real(r8) alpha, c1, c2
    real(r8) c(state%mesh%nlev+1)
    real(r8) A(state%mesh%nlev,3)
    real(r8) b(state%mesh%nlev)
    integer nlev, icol, k

    call calc_eddy_coef(state, dt)
    call mars_nasa_sfc_run(state)

    associate (mesh   => state%mesh   , &
               beta   => pbl_beta     , &
               dz     => state%dz     , & ! in
               u      => state%u      , & ! in
               v      => state%v      , & ! in
               rho    => state%rho    , & ! in
               ustar  => state%ustar  , & ! in
               tstar  => state%tstar  , & ! in
               cdh    => state%cdh    , & ! in
               km     => state%eddy_km, & ! in
               kh     => state%eddy_kh, & ! in
               taux   => state%taux   , & ! out
               tauy   => state%tauy   , & ! out
               hflx   => state%hflx   , & ! out
               rhouch => state%rhouch )   ! out
    nlev = mesh%nlev
    do icol = 1, mesh%ncol
      alpha = atan2(v(icol,nlev), u(icol,nlev))
      taux  (icol) =  rho(icol,nlev) * ustar(icol) * ustar(icol) * cos(alpha)
      tauy  (icol) =  rho(icol,nlev) * ustar(icol) * ustar(icol) * sin(alpha)
      hflx  (icol) = -rho(icol,nlev) * ustar(icol) * tstar(icol) * cpd
      rhouch(icol) =  rho(icol,nlev) * ustar(icol) * cdh  (icol) * cpd
      ! U
      do k = 2, mesh%nlev ! Interface levels excluding top and bottom
        c(k) = 2 * dt * km(icol,k) / (dz(icol,k-1) + dz(icol,k))
      end do
      do k = 1, mesh%nlev
        A(k,1) = -beta / dz(icol,k) * c(k  )
        A(k,3) = -beta / dz(icol,k) * c(k+1)
        A(k,2) = 1 - A(k,1) - A(k,3)
        c1 = (1 - beta) / dz(icol,k) * c(k  )
        c2 = (1 - beta) / dz(icol,k) * c(k+1)
        b(k) = c1 * u(icol,k-1) + (1 - c1 - c2) * u(icol,k) + c2 * u(icol,k+1)
      end do
    end do
    end associate

  end subroutine mars_nasa_pbl_run

  subroutine calc_eddy_coef(state, dt)

    type(mars_nasa_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    real(r8) rl, beta, dz, dptdz, dudz, dvdz, shear2, ri, km0, kh0, kmin
    integer icol, k

    associate (mesh         => state%mesh        , &
               z            => state%z           , & ! in
               z_lev        => state%z_lev       , & ! in
               pt           => state%pt          , & ! in
               u            => state%u           , & ! in
               v            => state%v           , & ! in
               saved_shear2 => state%saved_shear2, & ! inout
               km           => state%eddy_km     , & ! out
               kh           => state%eddy_kh     )   ! out
    do k = 2, mesh%nlev ! Interfaces excluding top and bottom.
      do icol = 1, mesh%ncol
        ! Calculate mixing length and beta (volume expansion coefficient).
        rl   = rl0 * ka * z_lev(icol,k) / (rl0 + ka * z_lev(icol,k))
        beta = (mesh%dlev(k-1) + mesh%dlev(k)) / (mesh%dlev(k-1) * pt(icol,k-1) + mesh%dlev(k) * pt(icol,k))
        ! Calculate gradient Richardson number.
        dz     = z(icol,k) - z(icol,k-1)
        dudz   = (u(icol,k) - u(icol,k-1)) / dz
        dvdz   = (v(icol,k) - v(icol,k-1)) / dz
        shear2 = dudz**2 + dvdz**2
        dptdz  = (pt(icol,k) - pt(icol,k-1)) / dz
        ! Smooth the wind shear used to calculate the gradient Richardson number.
        saved_shear2(icol,k) = saved_shear2(icol,k) - (saved_shear2(icol,k) - shear2) * dt / 1.0e4_r8
        ri = beta * g * dptdz / (saved_shear2(icol,k) + 1.0e-9_r8)
        ! Calculate neutral eddy coefficients.
        km0 = 0.393_r8 * rl**2 * sqrt(dudz / 0.153_r8)
        kh0 = 0.493_r8 * rl**2 * sqrt(dudz / 0.153_r8)
        ! Calculate eddy mixing coefficients.
        if (ri <= 0) then
          km(icol,k) = km0 * (1 - 15 * ri)**0.25_r8
          kh(icol,k) = kh0 * (1 - 15 * ri)**0.5_r8
        else
          km(icol,k) = km0 * (1 - ri / ric)
          kh(icol,k) = kh0 * (1 - ri / ric)
        end if
        ! Limit the coefficients.
        kmin = merge(0.1_r8, 0.001_r8, z(icol,k) < 300)
        km(icol,k) = max(km(icol,k), kmin)
        kh(icol,k) = max(kh(icol,k), kmin)
      end do
    end do
    end associate

  end subroutine calc_eddy_coef

end module mars_nasa_pbl_mod
