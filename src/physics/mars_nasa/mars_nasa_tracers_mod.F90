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

  implicit none

  private

  public mars_nasa_tracers_init
  public ntracers

  integer, public :: idx_m_dst = 0
  integer, public :: idx_n_dst = 0
  integer, public :: idx_m_vap = 0
  integer, public :: idx_m_cld = 0
  integer, public :: idx_n_cld = 0
  integer, public :: idx_m_ccn = 0

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

end module mars_nasa_tracers_mod
