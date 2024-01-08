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

  use flogger
  use string
  use process_mod

  implicit none

  character(1024) :: kcoef_file            = ''
  character(1024) :: dust_optics_file      = ''
  character(1024) :: cld_optics_file       = ''
  character(1024) :: albedo_file           = ''
  character(1024) :: thermal_inertia_file  = ''

  logical :: active_water   = .false.
  logical :: active_dust    = .false.

  integer :: nlev_soil = 0

  namelist /mars_nasa_control/ &
    kcoef_file               , &
    dust_optics_file         , &
    cld_optics_file          , &
    albedo_file              , &
    thermal_inertia_file     , &
    active_water             , &
    active_dust              , &
    nlev_soil

contains

  subroutine mars_nasa_parse_namelist(file_path, model_root)

    character(*), intent(in) :: file_path
    character(*), intent(in), optional :: model_root

    integer ierr
    logical is_exist

    open(10, file=file_path, status='old')
    read(10, nml=mars_nasa_control, iostat=ierr)
    if (ierr /= 0 .and. ierr /= -1) then
      call log_error('There is error in mars_nasa_control namelist in ' // trim(file_path) // '!', pid=proc%id)
    end if
    close(10)

    if (kcoef_file == '' .and. present(model_root)) then
      kcoef_file = abspath(trim(model_root) // '/data/mars/nasa_kcoef.nc')
    end if
    inquire(file=kcoef_file, exist=is_exist)
    if (.not. is_exist) then
      call log_error('kcoef_file ' // trim(kcoef_file) // ' does not exist!', pid=proc%id)
    end if

    if (dust_optics_file == '' .and. present(model_root)) then
      dust_optics_file = abspath(trim(model_root) // '/data/mars/nasa_dust_optics.nc')
    end if
    inquire(file=dust_optics_file, exist=is_exist)
    if (.not. is_exist) then
      call log_error('dust_optics_file ' // trim(dust_optics_file) // ' does not exist!', pid=proc%id)
    end if

    if (cld_optics_file == '' .and. present(model_root)) then
      cld_optics_file = abspath(trim(model_root) // '/data/mars/nasa_cld_optics.nc')
    end if
    inquire(file=cld_optics_file, exist=is_exist)
    if (.not. is_exist) then
      call log_error('cld_optics_file ' // trim(cld_optics_file) // ' does not exist!', pid=proc%id)
    end if

    if (albedo_file == '' .and. present(model_root)) then
      albedo_file = abspath(trim(model_root) // '/data/mars/mgs_albedo.nc')
    end if
    inquire(file=albedo_file, exist=is_exist)
    if (.not. is_exist) then
      call log_error('albedo_file ' // trim(albedo_file) // ' does not exist!', pid=proc%id)
    end if

    if (thermal_inertia_file == '' .and. present(model_root)) then
      thermal_inertia_file = abspath(trim(model_root) // '/data/mars/thermal_interia_2011.nc')
    end if
    inquire(file=thermal_inertia_file, exist=is_exist)
    if (.not. is_exist) then
      call log_error('thermal_inertia_file ' // trim(thermal_inertia_file) // ' does not exist!', pid=proc%id)
    end if

  end subroutine mars_nasa_parse_namelist

end module mars_nasa_namelist_mod
