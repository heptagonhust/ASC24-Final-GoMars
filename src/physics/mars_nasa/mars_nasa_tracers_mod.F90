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

module mars_nasa_tracers_mod

  ! Log-normal distribution
  !                1       ⎡ 1 ln(r / ro)²⎤
  !   n(r) dr = ------- exp⎢ - -----------⎥
  !             r 𝛔 √2𝜋    ⎣ 2     𝛔²     ⎦
  ! Median radius
  !        ⎡ 3 Mo  ⎤ 1/3
  !   ro = ⎢-------⎥     exp(-1.5 𝛔²)
  !        ⎣4𝜋 ⍴ No⎦ 
  ! Effective radius
  !
  !   reff = ro exp(2.5 𝛔²)

  use tracer_mod
  use mars_nasa_const_mod
  use mars_nasa_physics_types_mod

  implicit none

  private

  public mars_nasa_tracers_init

  integer, public :: idx_m_dst = 0
  integer, public :: idx_n_dst = 0
  integer, public :: idx_m_vap = 0
  integer, public :: idx_m_cld = 0
  integer, public :: idx_n_cld = 0
  integer, public :: idx_m_ccn = 0

  ! Dust particle density (kg m-3)
  real(r8), parameter :: rho_dst  = 2.5e3_r8
  ! Cloud ice particle density (kg m-3)
  real(r8), parameter :: rho_ice = 917.0_r8
  ! Standard deviation of dust particles
  real(r8), parameter :: dev2_dst = 0.63676_r8**2
  real(r8), parameter :: dev2_cld = 0.30870_r8**2

contains

  subroutine mars_nasa_tracers_init(dt_adv)

    real(r8), intent(in) :: dt_adv

    ! Water vapor is not size-dependent, so there is no water vapor number mixing ratio.
    ! The number of dust cores is assumed to equal to the number of ice particles.
    call tracer_add('mars', dt_adv, 'qm_dst', 'Dust mass mixing ratio'       , 'kg kg-1')
    idx_m_dst = ntracers
    call tracer_add('mars', dt_adv, 'qn_dst', 'Dust number mixing ratio'     , 'kg-1'   )
    idx_n_dst = ntracers
    call tracer_add('mars', dt_adv, 'qm_vap', 'Water vapor mass mixing ratio', 'kg kg-1')
    idx_m_vap = ntracers
    call tracer_add('mars', dt_adv, 'qm_cld', 'Ice cloud mass mixing ratio'  , 'kg kg-1')
    idx_m_cld = ntracers
    call tracer_add('mars', dt_adv, 'qn_cld', 'Ice cloud number mixing ratio', 'kg-1'   )
    idx_n_cld = ntracers
    call tracer_add('mars', dt_adv, 'qm_ccn', 'Dust core mass mixing ratio'  , 'kg kg-1')
    idx_m_ccn = ntracers

  end subroutine mars_nasa_tracers_init

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

end module mars_nasa_tracers_mod
