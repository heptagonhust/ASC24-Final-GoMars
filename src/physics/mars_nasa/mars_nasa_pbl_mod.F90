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

    real(r8) alpha
    integer icol

    call calc_eddy_coef(state, dt)
    call mars_nasa_sfc_run(state)

    associate (mesh   => state%mesh  , &
               u      => state%u     , & ! in
               v      => state%v     , & ! in
               rho    => state%rho   , & ! in
               ustar  => state%ustar , & ! in
               tstar  => state%tstar , & ! in
               cdh    => state%cdh   , & ! in
               taux   => state%taux  , & ! out
               tauy   => state%tauy  , & ! out
               hflx   => state%hflx  , & ! out
               rhouch => state%rhouch)   ! out
    do icol = 1, mesh%ncol
      alpha = atan2(v(icol,mesh%nlev), u(icol,mesh%nlev))
      taux  (icol) =  rho(icol,mesh%nlev) * ustar(icol) * ustar(icol) * cos(alpha)
      tauy  (icol) =  rho(icol,mesh%nlev) * ustar(icol) * ustar(icol) * sin(alpha)
      hflx  (icol) = -rho(icol,mesh%nlev) * ustar(icol) * tstar(icol) * cpd
      rhouch(icol) =  rho(icol,mesh%nlev) * ustar(icol) * cdh  (icol) * cpd
    end do
    end associate

  end subroutine mars_nasa_pbl_run

  subroutine calc_eddy_coef(state, dt)

    type(mars_nasa_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    real(r8) rl, beta, dz, dptdz, dudz, dvdz, shear2, ri, km0, kh0, kmin
    integer icol, k

    associate (mesh         => state%mesh        , &
               z            => state%z           , &
               z_lev        => state%z_lev       , &
               pt           => state%pt          , &
               u            => state%u           , &
               v            => state%v           , &
               saved_shear2 => state%saved_shear2, &
               km           => state%eddy_km     , &
               kh           => state%eddy_kh     )
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
