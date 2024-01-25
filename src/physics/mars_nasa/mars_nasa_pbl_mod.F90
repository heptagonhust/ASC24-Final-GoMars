module mars_nasa_pbl_mod

  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod

  implicit none

  private

  public mars_nasa_pbl_init
  public mars_nasa_pbl_final

contains

  subroutine mars_nasa_pbl_init()

  end subroutine mars_nasa_pbl_init

  subroutine mars_nasa_pbl_final()

  end subroutine mars_nasa_pbl_final

  subroutine mars_nasa_pbl_run()

  end subroutine mars_nasa_pbl_run

  subroutine calc_eddy_coef(state)

    type(mars_nasa_state_type), intent(inout) :: state

    real(r8) rl, dudz, ri, km0, kh0, kmin
    integer icol, k

    associate (mesh => state%mesh   , &
               z    => state%z      , &
               km   => state%eddy_km, &
               kh   => state%eddy_kh)
    do k = 1, mesh%nlev
      do icol = 1, mesh%ncol
        ! Calculate mixing length and beta.
        rl = rl0 * ka * z(icol,k) / (rl0 + ka * z(icol,k))
        ! Calculate Richardson number.

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
