module mars_nasa_physics_driver_mod

  use tracer_mod
  use mars_nasa_namelist_mod
  use mars_nasa_rad_mod

  implicit none

  private

  public mars_nasa_physics_driver_init
  public mars_nasa_physics_driver_final

contains

  subroutine mars_nasa_physics_driver_init(namelist_path, nlev)

    character(*), intent(in) :: namelist_path
    integer, intent(in) :: nlev

    call mars_nasa_parse_namelist(namelist_path)

    call mars_nasa_rad_init(nlev)

  end subroutine mars_nasa_physics_driver_init

  subroutine mars_nasa_physics_driver_final()

    call mars_nasa_rad_final()

  end subroutine mars_nasa_physics_driver_final

end module mars_nasa_physics_driver_mod
