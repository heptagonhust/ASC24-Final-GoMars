! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module const_mod

  use, intrinsic :: ieee_arithmetic
  use flogger
  use datetime, only: earth_day_seconds, mars_sol_seconds
  use gas_mod

  implicit none

#ifdef REAL_KIND
  integer, parameter :: r8 = REAL_KIND
#else
  integer, parameter :: r8 = 8
#endif

  real( 8), parameter :: pi     = atan(1.0d0) * 4.0d0
  real( 8), parameter :: pi2    = pi * 2
  real( 8), parameter :: pi05   = pi * 0.5d0
  real( 8), parameter :: deg    = 180.0d0 / pi
  real( 8), parameter :: rad    = pi / 180.0d0
  real(r8), parameter :: eps    = epsilon(1.0_r8)
  real(r8), parameter :: inf    = huge(1.0_r8)
  real(r8), parameter :: ka     = 0.4_r8             ! Karman constant

  real(r8)            :: omega        = 0 ! s-1
  real(r8)            :: radius       = 0 ! m
  real(r8)            :: periheli     = 0 ! Perihelion distance (million km)
  real(r8)            :: apheli       = 0 ! Aphelion distance (million km)
  real(r8)            :: eccen        = 0 ! Eccentricity
  real(r8)            :: obliq        = 0 ! Obliquity (deg)
  real(r8)            :: g            = 0 ! m2 s-2
  real(r8)            :: rd           = 0 ! J kg-1 K-1
  real(r8)            :: rv           = 0 ! J kg-1 K-1
  real(r8)            :: cpd          = 0 ! J kg-1 K-1
  real(r8)            :: cvd          = 0 ! J kg-1 K-1
  real(r8)            :: cpv          = 0 ! J kg-1 K-1
  real(r8)            :: cvv          = 0 ! J kg-1 K-1
  real(r8), parameter :: c_liq        = 4.1855d3 ! J kg-1 K-1
  real(r8), parameter :: c_ice        = 1.972d3  ! J kg-1 K-1
  real(r8)            :: lv           = 0 ! Latent heat of vaporization (J kg-1)
  real(r8)            :: rd_o_rv      = 0
  real(r8)            :: rv_o_rd      = 0
  real(r8)            :: rd_o_cpd     = 0
  real(r8)            :: rd_o_g       = 0
  real(r8)            :: cpd_o_cvd    = 0
  real(r8)            :: cvd_o_cpd    = 0
  real(r8)            :: lapse_rate   = 0 ! K m-1
  real(r8)            :: p0           = 0 ! Pa
  real(r8)            :: pk0          = 0 ! Pa

  integer, parameter :: inf_i4 = 10000000

  integer, parameter :: all_pass      = 0
  integer, parameter :: forward_pass  = 1
  integer, parameter :: backward_pass = 2
  integer, parameter :: nh_pass_1     = 3
  integer, parameter :: nh_pass_2     = 4
  integer            :: total_substeps = 0

  real(r8)           :: min_lon
  real(r8)           :: max_lon
  real(r8)           :: min_lat
  real(r8)           :: max_lat

contains

  subroutine const_init(planet)

    character(*), intent(in) :: planet

    rd  = major_gas%r
    cpd = major_gas%cp
    cvd = major_gas%cv

    select case (planet)
    case ('earth')
      omega      = 7.2921e-5_r8
      radius     = 6.37122d6
      periheli   = 147.1d0
      apheli     = 145.1d0
      obliq      = 23.4d0
      g          = 9.80616d0
      rv         = minor_gas%r
      cpv        = minor_gas%cp
      cvv        = minor_gas%cv
      lv         = minor_gas%l
      lapse_rate = 0.006d0
      rd_o_rv    = rd / rv
      rv_o_rd    = rv / rd
      p0         = 1.0d5
    case ('mars')
      omega      = 2 * pi / mars_sol_seconds
      radius     = 3.397200d6
      periheli   = 206.66d0
      apheli     = 249.22d0
      obliq      = 25.19d0
      g          = 3.72d0
      lv         = 2.84d6
      lapse_rate = 5.06d-3
      p0         = 6.1d2 ! FIXME: Should we use 6 hPa?
    case default
      call log_error('Invalid planet!')
    end select

    eccen     = (apheli - periheli) / (apheli + periheli)
    rd_o_g    = rd / g
    rd_o_cpd  = rd / cpd
    cpd_o_cvd = cpd / cvd
    cvd_o_cpd = cvd / cpd
    pk0       = p0**rd_o_cpd

  end subroutine const_init

end module const_mod
