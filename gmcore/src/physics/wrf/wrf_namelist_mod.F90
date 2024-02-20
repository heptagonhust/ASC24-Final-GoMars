module wrf_namelist_mod

  use const_mod

  implicit none

  integer :: num_soil_layers  = 0
  integer :: num_snow_layers  = 0
  integer :: urban_map_zd     = 0
  integer :: urban_map_zdf    = 0
  integer :: urban_map_zrd    = 0
  integer :: urban_map_zwd    = 0
  integer :: urban_map_zgrd   = 0
  integer :: urban_map_gd     = 0
  integer :: urban_map_bd     = 0
  integer :: urban_map_wd     = 0
  integer :: urban_map_gbd    = 0
  integer :: urban_map_fbd    = 0
  integer :: num_urban_ndm    = 0
  integer :: num_urban_hi     = 0

  namelist /wrf_control/  &
    num_soil_layers     , &
    num_snow_layers     , &
    urban_map_zd        , &
    urban_map_zdf       , &
    urban_map_zrd       , &
    urban_map_zwd       , &
    urban_map_zgrd      , &
    urban_map_gd        , &
    urban_map_bd        , &
    urban_map_wd        , &
    urban_map_gbd       , &
    urban_map_fbd       , &
    num_urban_ndm       , &
    num_urban_hi

end module wrf_namelist_mod
