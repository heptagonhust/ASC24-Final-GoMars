module mars_nasa_const_mod

  use const_mod
  use gas_mod, only: m_h2o, m_co2

  implicit none

  ! Surface albedo of ice at North Pole
  real(r8), parameter :: alb_ice_np     = 0.6_r8
  ! Surface albedo of ice at South Pole
  real(r8), parameter :: alb_ice_sp     = 0.5_r8
  ! Dust particle density (kg m-3)
  real(r8), parameter :: rho_dst        = 2.5e3_r8
  ! Cloud ice particle density (kg m-3)
  real(r8), parameter :: rho_ice        = 917.0_r8
  ! Standard deviation of dust particles
  real(r8), parameter :: dev2_dst       = 0.63676_r8**2
  real(r8), parameter :: dev2_cld       = 0.30870_r8**2
  ! Emisivity of bare ground in 15µm band, an absorption line of CO2
  real(r8), parameter :: emis_gnd_15um  = 1.0_r8
  ! A radiation code conversion factor
  real(r8), parameter :: cmk            = 3.51e+22_r8

  real(r8) :: ice_thresh_kgm2

end module mars_nasa_const_mod
