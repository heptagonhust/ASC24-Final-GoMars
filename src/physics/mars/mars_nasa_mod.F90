module mars_nasa_mod

  use const_mod

  implicit none

  real(r8), parameter :: ric_pbl  = 0.195_r8    ! Critical Richardson number for PBL
  real(r8), parameter :: factl    = 0.25_r8     ! Factor for defining soil layer thickness
  real(r8), parameter :: factm    = 1.2_r8      ! Factor for defining soil layer thickness
  real(r8), parameter :: skind    = 0.06_r8     ! Soil skin depth (m)
  real(r8), parameter :: npcwikg  = 0.0_r8      ! North polar cap water ice (kg)
  real(r8), parameter :: gidn     = 0.0545_r8   ! Ground ice depth in northern hemisphere (m)
  real(r8), parameter :: gids     = 0.0805_r8   ! Ground ice depth in southern hemisphere (m)

end module mars_nasa_mod