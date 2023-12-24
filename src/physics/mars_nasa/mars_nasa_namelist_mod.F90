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

module mars_nasa_namelist_mod

  implicit none

  character(256) kcoef_file
  character(256) dust_optics_file
  character(256) albedo_file
  character(256) thermal_inertia_file

  integer :: nlev_soil = 0

  namelist /mars_nasa_control/ &
    kcoef_file               , &
    dust_optics_file         , &
    albedo_file              , &
    thermal_inertia_file     , &
    nlev_soil

contains

  subroutine mars_nasa_parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=mars_nasa_control)
    close(10)

  end subroutine mars_nasa_parse_namelist

end module mars_nasa_namelist_mod
