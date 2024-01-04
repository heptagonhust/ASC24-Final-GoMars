module mars_nasa_const_mod

  use const_mod
  use gas_mod, only: m_h2o, m_co2

  implicit none

  ! A radiation code conversion factor
  real(r8), parameter :: cmk = 3.51e+22_r8

end module mars_nasa_const_mod
