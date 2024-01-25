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

module mars_nasa_physics_types_mod

  use fiona
  use process_mod
  use latlon_interp_mod
  use physics_types_mod
  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_spectra_mod
  use mars_nasa_rad_kcoef_mod
  use mars_nasa_tracers_mod

  implicit none

  private

  public mars_nasa_state_type
  public mars_nasa_tend_type
  public mars_nasa_static_type

  type, extends(physics_state_type) :: mars_nasa_state_type
    ! Surface latent heat flux (W m-2)
    real(r8), allocatable, dimension(    :    ) :: lhflx
    ! Atmosphere CO2 condensation (???)
    real(r8), allocatable, dimension(    :,:  ) :: atmcond
    ! Dust particle median radius
    real(r8), allocatable, dimension(    :,:  ) :: ro_dst
    ! Cloud particle median radius
    real(r8), allocatable, dimension(    :,:  ) :: ro_cld
    ! delta-Eddington optical thickness on the surface (1)
    real(r8), allocatable, dimension(:,:,:    ) :: detau
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_gas_vis
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_dst_vis
    real(r8), allocatable, dimension(:,:,:,:  ) :: tau_cld_vis
    ! Top or stratosphere temperature (K)
    real(r8), allocatable, dimension(      :  ) :: t_top
    real(r8), allocatable, dimension(      :  ) :: co2ice
    ! Boundary values of tracers (kg m-2???)
    real(r8), allocatable, dimension(    :,  :) :: q_sfc
    ! Intermediate variable ρ u* cp ch (???)
    real(r8), allocatable, dimension(    :    ) :: rhouch
    ! Eddy mixing coefficients in vertical (???)
    real(r8), allocatable, dimension(    :,:  ) :: eddy_km
    real(r8), allocatable, dimension(    :,:  ) :: eddy_kh
    ! Saved square of wind shear (s-2)
    real(r8), allocatable, dimension(    :,:  ) :: saved_shear2
    ! Heat exchange between atmosphere and surface (W m-2)
    real(r8), allocatable, dimension(    :    ) :: heat_sfc
    ! Soil conductivity (???)
    real(r8), allocatable, dimension(    :,:  ) :: soil_cond
    real(r8), allocatable, dimension(    :,:  ) :: soil_cond_lev
    ! Soil density (???)
    real(r8), allocatable, dimension(    :,:  ) :: soil_rho
    ! Soil specific heat capacity (???)
    real(r8), allocatable, dimension(    :,:  ) :: soil_cp
  contains
    procedure :: init  => mars_nasa_state_init
    procedure :: clear => mars_nasa_state_clear
    final mars_nasa_state_final
  end type mars_nasa_state_type

  type, extends(physics_tend_type) :: mars_nasa_tend_type
    real(r8), allocatable :: dpsdt(:)
  contains
    procedure :: init  => mars_nasa_tend_init
    procedure :: clear => mars_nasa_tend_clear
    final mars_nasa_tend_final
  end type mars_nasa_tend_type

  type, extends(physics_static_type) :: mars_nasa_static_type
    ! The fact that Mars is now a desert planet without oceans and lakes means
    ! that the thermal inertia of the surface is small.
    ! Surface thermal interia
    real(r8), allocatable, dimension(:  ) :: tin
    real(r8), allocatable, dimension(:,:) :: soil_tin
    ! Ground ice indicator for initializing soil model
    real(r8), allocatable, dimension(:  ) :: gnd_ice
  contains
    procedure :: init  => mars_nasa_static_init
    procedure :: clear => mars_nasa_static_clear
    procedure :: read  => mars_nasa_static_read
    final mars_nasa_static_final
  end type mars_nasa_static_type

contains

  subroutine mars_nasa_state_init(this, mesh)

    class(mars_nasa_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%lhflx        (                  mesh%ncol          ))
    allocate(this%atmcond      (                  mesh%ncol,mesh%nlev))
    allocate(this%ro_dst       (                  mesh%ncol,mesh%nlev))
    allocate(this%ro_cld       (                  mesh%ncol,mesh%nlev))
    allocate(this%detau        (spec_vis%n,ngauss,mesh%ncol          ))
    allocate(this%tau_gas_vis  (spec_vis%n,ngauss,mesh%ncol,mesh%nlev))
    allocate(this%tau_dst_vis  (spec_vis%n,ngauss,mesh%ncol,mesh%nlev))
    allocate(this%tau_cld_vis  (spec_vis%n,ngauss,mesh%ncol,mesh%nlev))
    allocate(this%t_top        (                  mesh%ncol          ))
    allocate(this%co2ice       (                  mesh%ncol          ))
    allocate(this%q_sfc        (                  mesh%ncol,         ntracers))
    allocate(this%rhouch       (                  mesh%ncol          ))
    allocate(this%eddy_km      (                  mesh%ncol,mesh%nlev+1))
    allocate(this%eddy_kh      (                  mesh%ncol,mesh%nlev+1))
    allocate(this%saved_shear2 (                  mesh%ncol,mesh%nlev+1))
    allocate(this%heat_sfc     (                  mesh%ncol          ))
    allocate(this%soil_cond    (                  mesh%ncol,nlev_soil))
    allocate(this%soil_cond_lev(                  mesh%ncol,nlev_soil+1))
    allocate(this%soil_rho     (                  mesh%ncol,nlev_soil))
    allocate(this%soil_cp      (                  mesh%ncol,nlev_soil))

    call this%physics_state_init(mesh)

  end subroutine mars_nasa_state_init

  subroutine mars_nasa_state_clear(this)

    class(mars_nasa_state_type), intent(inout) :: this

    if (allocated(this%lhflx        )) deallocate(this%lhflx        )
    if (allocated(this%atmcond      )) deallocate(this%atmcond      )
    if (allocated(this%ro_dst       )) deallocate(this%ro_dst       )
    if (allocated(this%ro_cld       )) deallocate(this%ro_cld       )
    if (allocated(this%detau        )) deallocate(this%detau        )
    if (allocated(this%tau_gas_vis  )) deallocate(this%tau_gas_vis  )
    if (allocated(this%tau_dst_vis  )) deallocate(this%tau_dst_vis  )
    if (allocated(this%tau_cld_vis  )) deallocate(this%tau_cld_vis  )
    if (allocated(this%t_top        )) deallocate(this%t_top        )
    if (allocated(this%co2ice       )) deallocate(this%co2ice       )
    if (allocated(this%q_sfc        )) deallocate(this%q_sfc        )
    if (allocated(this%rhouch       )) deallocate(this%rhouch       )
    if (allocated(this%eddy_km      )) deallocate(this%eddy_km      )
    if (allocated(this%eddy_kh      )) deallocate(this%eddy_kh      )
    if (allocated(this%saved_shear2 )) deallocate(this%saved_shear2 )
    if (allocated(this%heat_sfc     )) deallocate(this%heat_sfc     )
    if (allocated(this%soil_cond    )) deallocate(this%soil_cond    )
    if (allocated(this%soil_cond_lev)) deallocate(this%soil_cond_lev)
    if (allocated(this%soil_rho     )) deallocate(this%soil_rho     )
    if (allocated(this%soil_cp      )) deallocate(this%soil_cp      )

    call this%physics_state_clear()

  end subroutine mars_nasa_state_clear

  subroutine mars_nasa_state_final(this)

    type(mars_nasa_state_type), intent(inout) :: this

    call this%clear()

  end subroutine mars_nasa_state_final

  subroutine mars_nasa_tend_init(this, mesh)

    class(mars_nasa_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%dpsdt(mesh%ncol))

    call this%physics_tend_init(mesh)

  end subroutine mars_nasa_tend_init

  subroutine mars_nasa_tend_clear(this)

    class(mars_nasa_tend_type), intent(inout) :: this

    if (allocated(this%dpsdt)) deallocate(this%dpsdt)

    call this%physics_tend_clear()

  end subroutine mars_nasa_tend_clear

  subroutine mars_nasa_tend_final(this)

    type(mars_nasa_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine mars_nasa_tend_final

  subroutine mars_nasa_static_init(this, mesh)

    class(mars_nasa_static_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    call this%physics_static_init(mesh)

    allocate(this%tin     (mesh%ncol          ))
    allocate(this%soil_tin(mesh%ncol,nlev_soil))
    allocate(this%gnd_ice (mesh%ncol          ))

  end subroutine mars_nasa_static_init

  subroutine mars_nasa_static_clear(this)

    class(mars_nasa_static_type), intent(inout) :: this

    if (allocated(this%tin     )) deallocate(this%tin     )
    if (allocated(this%soil_tin)) deallocate(this%soil_tin)
    if (allocated(this%gnd_ice )) deallocate(this%gnd_ice )

    call this%physics_static_clear()

  end subroutine mars_nasa_static_clear

  subroutine mars_nasa_static_read(this, min_lon, max_lon, min_lat, max_lat, input_ngroup)

    class(mars_nasa_static_type), intent(inout) :: this
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat
    integer , intent(in) :: input_ngroup

    real(r8), allocatable :: lon(:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: array(:,:)

    ! Surface albedo
    call fiona_open_dataset('alb', file_path=albedo_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('alb', 'lon', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('alb', 'lat', span=[-90, 90])
    call fiona_start_input('alb')
    call fiona_input_range('alb', 'lon', lon, coord_range=[min_lon, max_lon]); lon = lon * rad
    call fiona_input_range('alb', 'lat', lat, coord_range=[min_lat, max_lat]); lat = lat * rad
    call fiona_input_range('alb', 'albedo', array, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_end_input('alb')
    call latlon_interp_bilinear_column(lon, lat, array, this%mesh%lon, this%mesh%lat, this%alb)
    deallocate(lon, lat, array)

    ! Surface thermal inertia
    call fiona_open_dataset('tin', file_path=thermal_inertia_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('tin', 'lon', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('tin', 'lat', span=[-90, 90])
    call fiona_start_input('tin')
    call fiona_input_range('tin', 'lon', lon, coord_range=[min_lon, max_lon]); lon = lon * rad
    call fiona_input_range('tin', 'lat', lat, coord_range=[min_lat, max_lat]); lat = lat * rad
    call fiona_input_range('tin', 'thin', array, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_end_input('tin')
    call latlon_interp_bilinear_column(lon, lat, array, this%mesh%lon, this%mesh%lat, this%tin)
    deallocate(lon, lat, array)

    ! Ground ice indicator
    call fiona_open_dataset('gnd_ice', file_path=gnd_ice_file, mpi_comm=proc%comm, ngroup=input_ngroup)
    call fiona_set_dim('gnd_ice', 'lon', span=[0, 360], cyclic=.true.)
    call fiona_set_dim('gnd_ice', 'lat', span=[-90, 90])
    call fiona_start_input('gnd_ice')
    call fiona_input_range('gnd_ice', 'lon', lon, coord_range=[min_lon, max_lon]); lon = lon * rad
    call fiona_input_range('gnd_ice', 'lat', lat, coord_range=[min_lat, max_lat]); lat = lat * rad
    call fiona_input_range('gnd_ice', 'gnd_ice', array, coord_range_1=[min_lon, max_lon], coord_range_2=[min_lat, max_lat])
    call fiona_end_input('gnd_ice')
    call latlon_interp_bilinear_column(lon, lat, array, this%mesh%lon, this%mesh%lat, this%gnd_ice)
    deallocate(lon, lat, array)

  end subroutine mars_nasa_static_read

  subroutine mars_nasa_static_final(this)

    type(mars_nasa_static_type), intent(inout) :: this

    call this%clear()

  end subroutine mars_nasa_static_final

end module mars_nasa_physics_types_mod
