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
  use mars_nasa_physics_types_mod
  use mars_nasa_physics_output_mod
  use mars_nasa_objects_mod
  use mars_orbit_mod
  use mars_nasa_solar_mod
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
      min_lat, max_lat, input_ngroup)

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
    integer , intent(in) :: input_ngroup

    call tracer_add('mars', dt_adv, 'qm_dst', 'Dust mass mixing ratio'           , 'kg kg-1'); idx_m_dst = ntracers
    call tracer_add('mars', dt_adv, 'qn_dst', 'Dust number mixing ratio'         , 'kg-1'   ); idx_n_dst = ntracers
    call tracer_add('mars', dt_adv, 'qm_vap', 'Water vapor mass mixing ratio'    , 'kg kg-1'); idx_m_vap = ntracers
    call tracer_add('mars', dt_adv, 'qm_cld', 'Cloud droplet mass mixing ratio'  , 'kg kg-1'); idx_m_cld = ntracers
    call tracer_add('mars', dt_adv, 'qn_cld', 'Cloud droplet number mixing ratio', 'kg-1'   ); idx_n_cld = ntracers
    call tracer_add('mars', dt_adv, 'qm_ccn', 'CCN mass mixing ratio'            , 'kg kg-1'); idx_m_ccn = ntracers

    call mars_nasa_physics_driver_final()
    call mars_nasa_objects_init(nblk, ncol, nlev, lon, lat, area)
    call mars_nasa_parse_namelist(namelist_path)
    call mars_nasa_read_static_data(min_lon, max_lon, min_lat, max_lat, input_ngroup)
    call mars_orbit_init()
    call mars_nasa_rad_init(nlev)

  end subroutine mars_nasa_physics_driver_init

  subroutine mars_nasa_physics_driver_final()

    call mars_nasa_rad_final()

    call mars_nasa_objects_final()

  end subroutine mars_nasa_physics_driver_final

  subroutine mars_nasa_physics_driver_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icol
    real(r8) ls

    ls = time%solar_longitude()

    call update_solar_decl_angle(ls)
    call update_solar_flux(ls)

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state)
      do icol = 1, objects(iblk)%mesh%ncol
        state%cosz(icol) = solar_cos_zenith_angle(mesh%lon(icol), mesh%lat(icol), time%hours_in_day())
      end do
      end associate
    end do

  end subroutine mars_nasa_physics_driver_run

  subroutine mars_nasa_physics_d2p()

  end subroutine mars_nasa_physics_d2p

  subroutine mars_nasa_physics_p2d()

  end subroutine mars_nasa_physics_p2d

end module mars_nasa_physics_driver_mod
