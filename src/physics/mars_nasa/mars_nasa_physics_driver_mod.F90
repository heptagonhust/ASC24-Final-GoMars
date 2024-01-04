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

module mars_nasa_physics_driver_mod

  use datetime
  use tracer_mod
  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_tracers_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_physics_output_mod
  use mars_nasa_objects_mod
  use mars_orbit_mod
  use mars_nasa_solar_mod
  use mars_nasa_rad_mod

  implicit none

  private

  public mars_nasa_init_stage2
  public mars_nasa_init_stage3
  public mars_nasa_final
  public mars_nasa_run
  public mars_nasa_d2p
  public mars_nasa_p2d
  public mars_nasa_add_output
  public mars_nasa_output
  public objects

  real(r8) dt

contains

  subroutine mars_nasa_init_stage2(namelist_path, mesh, dt_adv, dt_phys, &
    min_lon, max_lon, min_lat, max_lat, input_ngroup)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in), target :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat
    integer , intent(in) :: input_ngroup

    integer nlev

    nlev = mesh(1)%nlev

    call mars_nasa_final()
    call mars_nasa_parse_namelist(namelist_path)
    call mars_nasa_tracers_init(dt_adv)
    call mars_nasa_rad_init(nlev)
    call mars_nasa_objects_init(mesh)
    call mars_nasa_read_static_data(min_lon, max_lon, min_lat, max_lat, input_ngroup)
    call mars_orbit_init()

  end subroutine mars_nasa_init_stage2

  subroutine mars_nasa_init_stage3()

    integer iblk, icol

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do icol = 1, mesh%ncol
        state%t_top(icol) = state%t(icol,1)
        state%t_sfc(icol) = state%t(icol,mesh%nlev)
      end do
      end associate
    end do

  end subroutine mars_nasa_init_stage3

  subroutine mars_nasa_final()

    call mars_nasa_rad_final()

    call mars_nasa_objects_final()

  end subroutine mars_nasa_final

  subroutine mars_nasa_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icol
    real(r8) ls

    ls = time%solar_longitude()

    call update_solar_decl_angle(ls)
    call update_solar_flux(ls)

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
      do icol = 1, objects(iblk)%mesh%ncol
        state%cosz(icol) = solar_cos_zenith_angle(mesh%lon(icol), mesh%lat(icol), time%hours_in_day())
      end do
      call mars_nasa_rad_run(state, tend)
      end associate
    end do

  end subroutine mars_nasa_run

  subroutine mars_nasa_d2p()

  end subroutine mars_nasa_d2p

  subroutine mars_nasa_p2d()

  end subroutine mars_nasa_p2d

end module mars_nasa_physics_driver_mod
