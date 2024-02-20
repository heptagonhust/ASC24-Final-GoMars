module wrf_physics_driver_mod

  use wrf_namelist_mod
  use wrf_tracers_mod
  use wrf_physics_types_mod
  use wrf_objects_mod
  use lsm_driver_mod
  use pbl_driver_mod

  implicit none

  private

contains

  subroutine wrf_physics_init(namelist_path, mesh, dt_adv, dt_phys)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in) :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys

    call wrf_physics_final()
    call wrf_tracers_init(dt_adv)
    call wrf_objects_init(mesh)

  end subroutine wrf_physics_init

  subroutine wrf_physics_final()

    call wrf_objects_final()

  end subroutine wrf_physics_final

  subroutine wrf_physics_run()

  end subroutine wrf_physics_run

end module wrf_physics_driver_mod
