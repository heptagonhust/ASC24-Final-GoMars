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

  use fiona
  use tracer_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_objects_mod
  use mars_orbit_mod
  use mars_nasa_rad_mod

  implicit none

  private

  public mars_nasa_physics_driver_init
  public mars_nasa_physics_driver_final
  public mars_nasa_physics_driver_run
  public mars_nasa_physics_d2p
  public mars_nasa_physics_p2d
  public mars_nasa_physics_add_output
  public mars_nasa_physics_output
  public objects

  real(r8) dt

contains

  subroutine mars_nasa_physics_driver_init(namelist_path, nblk, &
      ncol, nlev, lon, lat, area, dt_adv, dt_phys, min_lon, max_lon, &
      min_lat, max_lat)

    character(*), intent(in) :: namelist_path
    integer , intent(in) :: nblk
    integer , intent(in) :: ncol(nblk)
    integer , intent(in) :: nlev
    real(r8), intent(in) :: lon (:,:) ! (ncol(iblk),nblk)
    real(r8), intent(in) :: lat (:,:) ! (ncol(iblk),nblk)
    real(r8), intent(in) :: area(:,:) ! (ncol(iblk),nblk)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat

    integer iblk

    call mars_nasa_physics_driver_final()
    call mars_nasa_objects_init(nblk, ncol, nlev, lon, lat, area)
    call mars_nasa_parse_namelist(namelist_path)
    call mars_nasa_read_static_data(min_lon, max_lon, min_lat, max_lat)
    call mars_orbit_init()
    call mars_nasa_rad_init(nlev)

  end subroutine mars_nasa_physics_driver_init

  subroutine mars_nasa_physics_driver_final()

    call mars_nasa_rad_final()

    call mars_nasa_objects_final()

  end subroutine mars_nasa_physics_driver_final

  subroutine mars_nasa_physics_driver_run()

  end subroutine mars_nasa_physics_driver_run

  subroutine mars_nasa_physics_d2p()

  end subroutine mars_nasa_physics_d2p

  subroutine mars_nasa_physics_p2d()

  end subroutine mars_nasa_physics_p2d

  subroutine mars_nasa_physics_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    call fiona_add_var(tag, 'alb', long_name='Surface albedo'         , units='', dim_names=['lon','lat'], dtype=dtype)
    call fiona_add_var(tag, 'tin', long_name='Surface thermal inertia', units='', dim_names=['lon','lat'], dtype=dtype)

  end subroutine mars_nasa_physics_add_output

  subroutine mars_nasa_physics_output(tag, iblk, start, count)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk
    integer, intent(in) :: start(3)
    integer, intent(in) :: count(3)

    associate (static => objects(iblk)%static)
    call fiona_output(tag, 'alb', reshape(static%sfc_alb, count(1:2)), start=start(1:2), count=count(1:2))
    call fiona_output(tag, 'tin', reshape(static%sfc_tin, count(1:2)), start=start(1:2), count=count(1:2))
    end associate

  end subroutine mars_nasa_physics_output

end module mars_nasa_physics_driver_mod
