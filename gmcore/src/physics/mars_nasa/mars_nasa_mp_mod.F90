module mars_nasa_mp_mod

  use mars_nasa_const_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_tracers_mod

  implicit none

  private

contains

  subroutine update_dust_params(state)

    type(mars_nasa_state_type), intent(inout) :: state

    real(r8), parameter :: f = 1.0_r8 / 3.0_r8
    real(r8), parameter :: c = 3.0_r8 / 4.0_r8 / pi
    real(r8) mo, no
    integer icol, k

    do k = 1, state%mesh%nlev
      do icol = 1, state%mesh%ncol
        mo = state%q(icol,k,idx_m_dst)
        no = state%q(icol,k,idx_n_dst)
        if (mo > 1.0e-8 .and. no > 1) then
          no = no + 1 ! FIXME: Why plus 1?
          state%ro_dst(icol,k) = (c * mo / no / rho_dst)**f * exp(-1.5_r8 * dev2_dst)
        end if
      end do
    end do

  end subroutine update_dust_params

  subroutine update_cloud_params(state)

    type(mars_nasa_state_type), intent(inout) :: state

    real(r8), parameter :: f = 1.0_r8 / 3.0_r8
    real(r8), parameter :: c = 3.0_r8 / 4.0_r8 / pi
    real(r8) mo, no, rho_cld
    integer icol, k

    do k = 1, state%mesh%nlev
      do icol = 1, state%mesh%ncol
        mo = state%q(icol,k,idx_m_cld) + state%q(icol,k,idx_m_ccn)
        no = state%q(icol,k,idx_n_cld)
        if (mo > 1.0e-7 .and. no > 1) then
          rho_cld = state%q(icol,k,idx_m_cld) / mo * rho_ice + &
                    state%q(icol,k,idx_m_ccn) / mo * rho_dst
          state%ro_cld(icol,k) = (c * mo / no / rho_cld)**f * exp(-1.5_r8 * dev2_cld)
        end if
      end do
    end do

  end subroutine update_cloud_params

end module mars_nasa_mp_mod