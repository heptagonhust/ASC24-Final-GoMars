module mars_nasa_physics_driver_mod

  use namelist_mod
  use tracer_mod
  use mars_nasa_rad_mod

  implicit none

  private

  public mars_nasa_physics_driver_init
  public mars_nasa_physics_driver_final

  character(256) kcoef_file_path

  namelist /mars_nasa_control/ &
    kcoef_file_path

contains

  subroutine mars_nasa_physics_driver_init(namelist_path, nlev)

    character(*), intent(in) :: namelist_path
    integer, intent(in) :: nlev

    open(10, file=namelist_path, status='old')
    read(10, nml=mars_nasa_control)
    close(10)

    call tracer_add('mars', dt_adv, 'qd_m', 'kg kg-1')
    call tracer_add('mars', dt_adv, 'qd_n', 'kg-1')
    call tracer_add('mars', dt_adv, 'qc_m', 'kg kg-1')
    call tracer_add('mars', dt_adv, 'qc_n', 'kg-1')
    call tracer_add('mars', dt_adv, 'qdc' , 'kg kg-1')
    call tracer_add('mars', dt_adv, 'qv_m', 'kg kg-1')

    call mars_nasa_rad_init(nlev, kcoef_file_path)

  end subroutine mars_nasa_physics_driver_init

  subroutine mars_nasa_physics_driver_final()

    call mars_nasa_rad_final()

  end subroutine mars_nasa_physics_driver_final

end module mars_nasa_physics_driver_mod
