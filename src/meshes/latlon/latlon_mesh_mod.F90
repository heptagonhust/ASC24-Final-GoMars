module latlon_mesh_mod

  use flogger
  use const_mod, only: pi, pi2, pi05, radius, inf, deg
  use sphere_geometry_mod

  implicit none

  private

  public latlon_mesh_type
  public global_mesh

  type latlon_mesh_type
    ! For nesting
    integer :: id = 0
    type(latlon_mesh_type), pointer :: parent => null()
    integer lon_hw              ! Halo width along longitude
    integer lat_hw              ! Halo width along latitude
    integer full_nlon
    integer half_nlon
    integer full_nlat
    integer half_nlat
    integer full_nlev
    integer half_nlev
    integer full_ids, full_ide
    integer full_jds, full_jde
    integer full_jds_no_pole
    integer full_jde_no_pole
    integer full_kds, full_kde
    integer half_ids, half_ide
    integer half_jds, half_jde
    integer half_kds, half_kde
    integer full_ims, full_ime
    integer half_ims, half_ime
    integer full_jms, full_jme
    integer half_jms, half_jme
    integer full_kms, full_kme
    integer half_kms, half_kme
    real(8) start_lon
    real(8) end_lon
    real(8) start_lat
    real(8) end_lat
    real(8) dlon
    real(8), allocatable, dimension(:  ) :: dlat
    real(8), allocatable, dimension(:  ) :: full_dlev
    real(8), allocatable, dimension(:  ) :: half_dlev
    real(8), allocatable, dimension(:  ) :: half_dlev_upper
    real(8), allocatable, dimension(:  ) :: half_dlev_lower
    real(8) total_area
    real(8), allocatable, dimension(:  ) :: full_lon
    real(8), allocatable, dimension(:  ) :: half_lon
    real(8), allocatable, dimension(:  ) :: full_lat
    real(8), allocatable, dimension(:  ) :: half_lat
    real(8), allocatable, dimension(:  ) :: full_lev
    real(8), allocatable, dimension(:  ) :: half_lev
    real(8), allocatable, dimension(:  ) :: full_cos_lon
    real(8), allocatable, dimension(:  ) :: half_cos_lon
    real(8), allocatable, dimension(:  ) :: full_sin_lon
    real(8), allocatable, dimension(:  ) :: half_sin_lon
    real(8), allocatable, dimension(:  ) :: full_cos_lat
    real(8), allocatable, dimension(:  ) :: half_cos_lat
    real(8), allocatable, dimension(:  ) :: full_sin_lat
    real(8), allocatable, dimension(:  ) :: half_sin_lat
    ! For output
    real(8), allocatable, dimension(:  ) :: full_lon_deg
    real(8), allocatable, dimension(:  ) :: half_lon_deg
    real(8), allocatable, dimension(:  ) :: full_lat_deg
    real(8), allocatable, dimension(:  ) :: half_lat_deg
    ! Area for weighting
    real(8), allocatable, dimension(:  ) :: area_cell
    real(8), allocatable, dimension(:  ) :: area_lon
    real(8), allocatable, dimension(:  ) :: area_lon_west
    real(8), allocatable, dimension(:  ) :: area_lon_east
    real(8), allocatable, dimension(:  ) :: area_lon_north
    real(8), allocatable, dimension(:  ) :: area_lon_south
    real(8), allocatable, dimension(:  ) :: area_lat
    real(8), allocatable, dimension(:  ) :: area_lat_west
    real(8), allocatable, dimension(:  ) :: area_lat_east
    real(8), allocatable, dimension(:  ) :: area_lat_north
    real(8), allocatable, dimension(:  ) :: area_lat_south
    real(8), allocatable, dimension(:  ) :: area_vtx
    real(8), allocatable, dimension(:,:) :: area_subcell
    ! Edge length
    real(8), allocatable, dimension(:  ) :: de_lon
    real(8), allocatable, dimension(:  ) :: de_lat
    real(8), allocatable, dimension(:  ) :: le_lat
    real(8), allocatable, dimension(:  ) :: le_lon
  contains
    procedure :: init_global                  => latlon_mesh_init_global
    procedure :: init_from_parent             => latlon_mesh_init_from_parent
    procedure :: reinit                       => latlon_mesh_reinit
    procedure :: common_init                  => latlon_mesh_common_init
    procedure :: has_south_pole               => latlon_mesh_has_south_pole
    procedure :: has_north_pole               => latlon_mesh_has_north_pole
    procedure :: is_south_pole                => latlon_mesh_is_south_pole
    procedure :: is_north_pole                => latlon_mesh_is_north_pole
    procedure :: is_pole                      => latlon_mesh_is_pole
    procedure :: is_full_lat_next_to_pole     => latlon_mesh_is_full_lat_next_to_pole
    procedure :: is_half_lat_next_to_pole     => latlon_mesh_is_half_lat_next_to_pole
    procedure :: is_inside_with_halo_full_lat => latlon_mesh_is_inside_with_halo_full_lat
    procedure :: is_inside_with_halo_half_lat => latlon_mesh_is_inside_with_halo_half_lat
    procedure :: is_inside_pole_full_lat      => latlon_mesh_is_inside_pole_full_lat
    procedure :: is_inside_pole_half_lat      => latlon_mesh_is_inside_pole_half_lat
    procedure :: is_outside_pole_full_lat     => latlon_mesh_is_outside_pole_full_lat
    procedure :: is_outside_pole_half_lat     => latlon_mesh_is_outside_pole_half_lat
    procedure :: clear                        => latlon_mesh_clear
    final :: latlon_mesh_final
  end type latlon_mesh_type

  type(latlon_mesh_type), target :: global_mesh

contains

  subroutine latlon_mesh_init_global(this, nlon, nlat, nlev, id, lon_hw, lat_hw, coarse_pole_decay, coarse_pole_mul, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    integer, intent(in)           :: nlon
    integer, intent(in)           :: nlat
    integer, intent(in), optional :: nlev
    integer, intent(in), optional :: id
    integer, intent(in), optional :: lon_hw
    integer, intent(in), optional :: lat_hw
    real(8), intent(in), optional :: coarse_pole_decay
    real(8), intent(in), optional :: coarse_pole_mul
    logical, intent(in), optional :: keep_lev

    integer nlev_opt, id_opt, lon_hw_opt, lat_hw_opt
    real(8) coarse_pole_decay_opt, coarse_pole_mul_opt
    logical keep_lev_opt
    real(8) dlat0
    real(16) x(3), y(3), z(3)
    integer i, j, ierr

    nlev_opt              =  1; if (present(nlev             )) nlev_opt              = nlev
    id_opt                = -1; if (present(id               )) id_opt                = id
    lon_hw_opt            =  3; if (present(lon_hw           )) lon_hw_opt            = lon_hw
    lat_hw_opt            =  3; if (present(lat_hw           )) lat_hw_opt            = lat_hw
    coarse_pole_decay_opt =  0; if (present(coarse_pole_decay)) coarse_pole_decay_opt = coarse_pole_decay
    coarse_pole_mul_opt   =  0; if (present(coarse_pole_mul  )) coarse_pole_mul_opt   = coarse_pole_mul

    call this%clear(keep_lev)

    this%full_nlon = nlon
    this%half_nlon = nlon
    this%full_ids  = 1
    this%full_ide  = this%full_nlon
    this%half_ids  = 1
    this%half_ide  = this%half_nlon
    this%full_nlat = nlat
    this%half_nlat = nlat - 1
    this%full_jds  = 1
    this%full_jde  = this%full_nlat
    this%half_jds  = 1
    this%half_jde  = this%half_nlat
    this%full_nlev = nlev_opt
    this%half_nlev = this%full_nlev + 1
    this%full_kds  = 1
    this%full_kde  = this%full_nlev
    this%half_kds  = 1
    this%half_kde  = this%half_nlev

    this%id        = id_opt
    this%lon_hw    = lon_hw_opt
    this%lat_hw    = lat_hw_opt
    this%start_lon = 0
    this%end_lon   =  pi2
    this%start_lat = -pi05
    this%end_lat   =  pi05

    call this%common_init()

    this%dlon = (this%end_lon - this%start_lon) / this%full_nlon
    do i = this%full_ims, this%full_ime
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5d0 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

    ! Set initial guess latitudes of full merdional grids.
    dlat0 = (this%end_lat - this%start_lat) / this%half_nlat
    do j = 1, this%half_nlat
      this%half_lat(j) = this%start_lat + (j - 0.5d0) * dlat0
      if (abs(this%half_lat(j)) < 1.0e-12) this%half_lat(j) = 0
    end do

    if (coarse_pole_mul_opt /= 0) then
      ! Calculate real dlat which is large at polar region.
      dlat0 = this%dlon
      do j = 1, this%half_nlat
        this%dlat(j) = dlat0 * (1 + (coarse_pole_mul_opt - 1) * exp(-coarse_pole_decay_opt * (abs(this%half_lat(j)) - pi05)**2))
      end do
      this%dlat(1:this%half_nlat) = this%dlat(1:this%half_nlat) / sum(this%dlat(1:this%half_nlat)) * pi
    else
      this%dlat(1:this%half_nlat) = dlat0
    end if

    ! Set latitudes of full merdional grids.
    this%full_lat(1) = this%start_lat
    this%full_lat_deg(1) = this%start_lat * deg
    do j = 2, this%full_nlat - 1
      this%full_lat(j) = this%full_lat(j-1) + this%dlat(j-1)
      if (abs(this%full_lat(j)) < 1.0e-12) this%full_lat(j) = 0
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
    this%full_lat(this%full_nlat) = this%end_lat
    this%full_lat_deg(this%full_nlat) = this%end_lat * deg

    ! Set latitudes of half merdional grids.
    do j = 1, this%half_nlat
      if (this%full_lat(j) == pi05) cycle
      this%half_lat(j) = this%full_lat(j) + 0.5d0 * this%dlat(j)
      if (abs(this%half_lat(j)) < 1.0e-12) this%half_lat(j) = 0
      this%half_lat_deg(j) = this%half_lat(j) * deg
    end do

    ! Ensure the grids are equatorial symmetry.
    do j = 1, this%full_nlat
      if (this%full_lat(j) > 0) then
        this%full_lat(j) = -this%full_lat(this%full_nlat-j+1)
        this%full_lat_deg(j) = -this%full_lat_deg(this%full_nlat-j+1)
      end if
    end do
    do j = 1, this%half_nlat
      if (this%half_lat(j) > 0) then
        this%half_lat(j) = -this%half_lat(this%half_nlat-j+1)
        this%half_lat_deg(j) = -this%half_lat_deg(this%half_nlat-j+1)
      end if
    end do

    do i = this%full_ims, this%full_ime
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do

    do i = this%half_ims, this%half_ime
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do

    do j = this%half_jms, this%half_jme
      if (this%half_lat(j) >= -pi05 .and. this%half_lat(j) <= pi05) then
        this%half_cos_lat(j) = cos(this%half_lat(j))
        this%half_sin_lat(j) = sin(this%half_lat(j))
      end if
    end do

    do j = this%full_jms, this%full_jme
      if (this%full_lat(j) >= -pi05 .and. this%full_lat(j) <= pi05) then
        this%full_cos_lat(j) = cos(this%full_lat(j))
        this%full_sin_lat(j) = sin(this%full_lat(j))
      end if
    end do

    ! Ensure the values of cos_lat and sin_lat are expected at the Poles.
    this%full_cos_lat(this%full_jds) =  0
    this%full_sin_lat(this%full_jds) = -1
    this%full_cos_lat(this%full_jde) =  0
    this%full_sin_lat(this%full_jde) =  1

    do j = this%full_jds, this%full_jde
      if (this%is_south_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
      else if (this%is_north_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (1.0d0 - this%half_sin_lat(j-1))
      else
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        !
        !           1,j
        !           /|
        !          / |
        !         /  |
        !        /   |
        !    1,j \   |
        !         \  |
        !          \ |
        !           \|
        !          1,j-1
        !
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j-1), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%area_lon_west(j) = spherical_area(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
        this%area_lon_east(j) = this%area_lon_west(j)
        this%area_lon(j) = this%area_lon_west(j) + this%area_lon_east(j)
        !
        !          1,j
        !           /\
        !          /  \
        !         /    \
        !        /______\
        !    1,j          2,j
        !
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_north(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
        !
        !    1,j          2,j
        !        --------
        !        \      /
        !         \    /
        !          \  /
        !           \/
        !         1,j-1
        !
        call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j-1), x(1), y(1), z(1))
        call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
        call lonlat2xyz(radius, this%full_lon(1), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_south(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
        if (ierr /= 0) then
          call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
        end if
      end if
    end do

    do j = this%half_jds, this%half_jde
      !
      !          2,j+1
      !           /|
      !          / |
      !         /  |
      !        /   |
      !    1,j \   |
      !         \  |
      !          \ |
      !           \|
      !           2,j
      !
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j+1), x(3), y(3), z(3))
      this%area_lat_west(j) = spherical_area(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      this%area_lat_east(j) = this%area_lat_west(j)
      !
      !         2,j+1
      !           /\
      !          /  \
      !         /    \
      !        /______\
      !    1,j          2,j
      !
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j+1), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
      this%area_lat_north(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      !
      !    1,j          2,j
      !        --------
      !        \      /
      !         \    /
      !          \  /
      !           \/
      !          2,j
      !
      call lonlat2xyz(radius, this%full_lon(2), this%full_lat(j), x(1), y(1), z(1))
      call lonlat2xyz(radius, this%half_lon(2), this%half_lat(j), x(2), y(2), z(2))
      call lonlat2xyz(radius, this%half_lon(1), this%half_lat(j), x(3), y(3), z(3))
      this%area_lat_south(j) = spherical_area_with_last_small_arc(radius, x, y, z, ierr)
      if (ierr /= 0) then
        call log_error(sphere_geometry_error_message(ierr), __FILE__, __LINE__)
      end if
      ! Reset up or down area to polar sector area.
      if (this%is_south_pole(j)) then
        this%area_lat_south(j) = this%area_cell(j)
      else if (this%is_north_pole(j+1)) then
        this%area_lat_north(j) = this%area_cell(j+1)
      end if
      this%area_lat(j) = this%area_lat_north(j) + this%area_lat_south(j)
    end do

    do j = this%half_jds, this%half_jde
      if (this%is_south_pole(j)) then
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_south(j+1)
      else if (this%is_north_pole(j+1)) then
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_north(j)
      else
        this%area_vtx(j) = this%area_lat_west(j) + this%area_lat_east(j) + this%area_lon_south(j+1) + this%area_lon_north(j)
      end if
    end do

    do j = this%full_jds_no_pole, this%full_jde_no_pole
      this%de_lon(j) = radius * this%full_cos_lat(j) * this%dlon
      this%le_lon(j) = 2.0d0 * this%area_lon(j) / this%de_lon(j)
    end do
    if (this%has_south_pole()) then
      this%le_lon(this%full_jds) = 0
      this%de_lon(this%full_jds) = 0
    end if
    if (this%has_north_pole()) then
      this%le_lon(this%full_jde) = 0
      this%de_lon(this%full_jde) = 0
    end if

    do j = this%half_jds, this%half_jde
      this%le_lat(j) = radius * this%half_cos_lat(j) * this%dlon
      this%de_lat(j) = 2.0d0 * this%area_lat(j) / this%le_lat(j)
    end do

  end subroutine latlon_mesh_init_global

  subroutine latlon_mesh_init_from_parent(this, parent, id, ids, ide, jds, jde, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    class(latlon_mesh_type), intent(in), target :: parent
    integer, intent(in) :: id
    integer, intent(in) :: ids
    integer, intent(in) :: ide
    integer, intent(in) :: jds
    integer, intent(in) :: jde
    logical, intent(in), optional :: keep_lev

    integer i, j

    call this%clear(keep_lev)

    this%parent => parent

    this%full_nlon = ide - ids + 1
    this%half_nlon = this%full_nlon
    this%full_ids  = ids
    this%full_ide  = ide
    this%half_ids  = ids
    this%half_ide  = ide
    this%full_nlat = jde - jds + 1
    this%full_jds  = jds
    this%full_jde  = jde
    this%half_jds  = jds
    this%half_jde  = merge(jde - 1, jde, this%has_north_pole())
    this%half_nlat = this%half_jde - this%half_jds + 1

    this%full_nlev = parent%full_nlev
    this%half_nlev = parent%half_nlev
    this%full_kds  = parent%full_kds
    this%full_kde  = parent%full_kde
    this%half_kds  = parent%half_kds
    this%half_kde  = parent%half_kde

    this%id        = id
    this%lon_hw    = parent%lon_hw
    this%lat_hw    = parent%lat_hw
    this%start_lon = parent%full_lon(ids)
    this%end_lon   = parent%full_lon(ide) + parent%dlon
    this%start_lat = merge(parent%half_lat(jds-1), -pi05, .not. this%has_south_pole())
    this%end_lat   = merge(parent%half_lat(jde  ),  pi05, .not. this%has_north_pole())

    call this%common_init()

    this%full_dlev = parent%full_dlev
    this%half_dlev = parent%half_dlev
    this%half_dlev_upper = parent%half_dlev_upper
    this%half_dlev_lower = parent%half_dlev_lower

    this%dlon = parent%dlon
    do i = this%full_ims, this%full_ime
      this%full_lon(i)     = parent%full_lon(i)
      this%half_lon(i)     = parent%half_lon(i)
      this%full_lon_deg(i) = parent%full_lon_deg(i)
      this%half_lon_deg(i) = parent%half_lon_deg(i)
      this%full_sin_lon(i) = parent%full_sin_lon(i)
      this%half_sin_lon(i) = parent%half_sin_lon(i)
      this%full_cos_lon(i) = parent%full_cos_lon(i)
      this%half_cos_lon(i) = parent%half_cos_lon(i)
    end do

    this%dlat = parent%dlat(lbound(this%dlat, 1):ubound(this%dlat, 1))
    do j = this%full_jms, this%full_jme
      this%full_lat(j)       = parent%full_lat(j)
      this%full_lat_deg(j)   = parent%full_lat_deg(j)
      this%full_sin_lat(j)   = parent%full_sin_lat(j)
      this%full_cos_lat(j)   = parent%full_cos_lat(j)
      this%area_cell(j)      = parent%area_cell(j)
      this%area_subcell(:,j) = parent%area_subcell(:,j)
      this%area_lon_west(j)  = parent%area_lon_west(j)
      this%area_lon_east(j)  = parent%area_lon_east(j)
      this%area_lon_north(j) = parent%area_lon_north(j)
      this%area_lon_south(j) = parent%area_lon_south(j)
      this%area_lon(j)       = parent%area_lon(j)
      this%le_lon(j)         = parent%le_lon(j)
      this%de_lon(j)         = parent%de_lon(j)
    end do
    do j = this%half_jms, this%half_jme
      this%half_lat(j)       = parent%half_lat(j)
      this%half_lat_deg(j)   = parent%half_lat_deg(j)
      this%half_sin_lat(j)   = parent%half_sin_lat(j)
      this%half_cos_lat(j)   = parent%half_cos_lat(j)
      this%area_vtx(j)       = parent%area_vtx(j)
      this%area_lat_west(j)  = parent%area_lat_west(j)
      this%area_lat_east(j)  = parent%area_lat_east(j)
      this%area_lat_north(j) = parent%area_lat_north(j)
      this%area_lat_south(j) = parent%area_lat_south(j)
      this%area_lat(j)       = parent%area_lat(j)
      this%le_lat(j)         = parent%le_lat(j)
      this%de_lat(j)         = parent%de_lat(j)
    end do

    this%full_lev = parent%full_lev
    this%half_lev = parent%half_lev

  end subroutine latlon_mesh_init_from_parent

  subroutine latlon_mesh_reinit(this, lon_hw)

    class(latlon_mesh_type), intent(inout) :: this
    integer, intent(in), optional :: lon_hw

    integer nlon, nlat, nlev, lat_hw

    nlon = global_mesh%full_nlon
    nlat = global_mesh%full_nlat
    nlev = global_mesh%full_nlev
    lat_hw = global_mesh%lat_hw

    if (associated(this%parent)) then
      call this%init_from_parent(this%parent, this%id, this%full_ids, this%full_ide, this%full_jds, this%full_jde, keep_lev=.true.)
    else if (present(lon_hw)) then
      call this%init_global(nlon, nlat, nlev, 0, lon_hw, lat_hw, keep_lev=.true.)
    else
      call log_error('Logical error!', __FILE__, __LINE__)
    end if

  end subroutine latlon_mesh_reinit

  subroutine latlon_mesh_common_init(this)

    class(latlon_mesh_type), intent(inout) :: this

    this%total_area = radius**2 * (this%end_lon - this%start_lon) * (sin(this%end_lat) - sin(this%start_lat))

    this%full_jds_no_pole = merge(this%full_jds + 1, this%full_jds, this%has_south_pole())
    this%full_jde_no_pole = merge(this%full_jde - 1, this%full_jde, this%has_north_pole())
    this%half_jds = this%half_jds
    this%half_jde = this%half_jde

    ! Use maximum lon_hw in this process and its south and north neighbors.
    this%full_ims = this%full_ids - this%lon_hw
    this%full_ime = this%full_ide + this%lon_hw
    this%full_jms = this%full_jds - this%lat_hw
    this%full_jme = this%full_jde + this%lat_hw
    this%half_ims = this%half_ids - this%lon_hw
    this%half_ime = this%half_ide + this%lon_hw
    this%half_jms = this%half_jds - this%lat_hw
    this%half_jme = this%half_jde + this%lat_hw
    this%full_kms = this%full_kds - 3
    this%full_kme = this%full_kde + 3
    this%half_kms = this%half_kds
    this%half_kme = this%half_kde

    if (.not. allocated(this%full_lev)) then
      allocate(this%full_dlev        (this%full_kms:this%full_kme)); this%full_dlev           = 0
      allocate(this%half_dlev        (this%half_kms:this%half_kme)); this%half_dlev           = 0
      allocate(this%half_dlev_upper  (this%half_kms:this%half_kme)); this%half_dlev_upper     = 0
      allocate(this%half_dlev_lower  (this%half_kms:this%half_kme)); this%half_dlev_lower     = 0
      allocate(this%full_lev         (this%full_kms:this%full_kme)); this%full_lev            = inf
      allocate(this%half_lev         (this%half_kms:this%half_kme)); this%half_lev            = inf
    end if

    allocate(this%dlat               (this%half_jms:this%half_jme)); this%dlat                = 0
    allocate(this%full_lon           (this%full_ims:this%full_ime)); this%full_lon            = inf
    allocate(this%half_lon           (this%half_ims:this%half_ime)); this%half_lon            = inf
    allocate(this%full_lat           (this%full_jms:this%full_jme)); this%full_lat            = inf
    allocate(this%half_lat           (this%half_jms:this%half_jme)); this%half_lat            = inf
    allocate(this%full_cos_lon       (this%full_ims:this%full_ime)); this%full_cos_lon        = inf
    allocate(this%half_cos_lon       (this%half_ims:this%half_ime)); this%half_cos_lon        = inf
    allocate(this%full_sin_lon       (this%full_ims:this%full_ime)); this%full_sin_lon        = inf
    allocate(this%half_sin_lon       (this%half_ims:this%half_ime)); this%half_sin_lon        = inf
    allocate(this%full_cos_lat       (this%full_jms:this%full_jme)); this%full_cos_lat        = inf
    allocate(this%half_cos_lat       (this%half_jms:this%half_jme)); this%half_cos_lat        = inf
    allocate(this%full_sin_lat       (this%full_jms:this%full_jme)); this%full_sin_lat        = inf
    allocate(this%half_sin_lat       (this%half_jms:this%half_jme)); this%half_sin_lat        = inf
    allocate(this%full_lon_deg       (this%full_ims:this%full_ime)); this%full_lon_deg        = inf
    allocate(this%half_lon_deg       (this%half_ims:this%half_ime)); this%half_lon_deg        = inf
    allocate(this%full_lat_deg       (this%full_jms:this%full_jme)); this%full_lat_deg        = inf
    allocate(this%half_lat_deg       (this%half_jms:this%half_jme)); this%half_lat_deg        = inf
    allocate(this%area_cell          (this%full_jms:this%full_jme)); this%area_cell           = 0
    allocate(this%area_lon           (this%full_jms:this%full_jme)); this%area_lon            = 0
    allocate(this%area_lon_west      (this%full_jms:this%full_jme)); this%area_lon_west       = 0
    allocate(this%area_lon_east      (this%full_jms:this%full_jme)); this%area_lon_east       = 0
    allocate(this%area_lon_north     (this%full_jms:this%full_jme)); this%area_lon_north      = 0
    allocate(this%area_lon_south     (this%full_jms:this%full_jme)); this%area_lon_south      = 0
    allocate(this%area_lat           (this%half_jms:this%half_jme)); this%area_lat            = 0
    allocate(this%area_lat_west      (this%half_jms:this%half_jme)); this%area_lat_west       = 0
    allocate(this%area_lat_east      (this%half_jms:this%half_jme)); this%area_lat_east       = 0
    allocate(this%area_lat_north     (this%half_jms:this%half_jme)); this%area_lat_north      = 0
    allocate(this%area_lat_south     (this%half_jms:this%half_jme)); this%area_lat_south      = 0
    allocate(this%area_vtx           (this%half_jms:this%half_jme)); this%area_vtx            = 0
    allocate(this%area_subcell     (2,this%full_jms:this%full_jme)); this%area_subcell        = 0
    allocate(this%de_lon             (this%full_jms:this%full_jme)); this%de_lon              = 0
    allocate(this%de_lat             (this%half_jms:this%half_jme)); this%de_lat              = 0
    allocate(this%le_lat             (this%half_jms:this%half_jme)); this%le_lat              = 0
    allocate(this%le_lon             (this%full_jms:this%full_jme)); this%le_lon              = 0

  end subroutine latlon_mesh_common_init

  pure logical function latlon_mesh_has_south_pole(this) result(res)

    class(latlon_mesh_type), intent(in) :: this

    res = this%full_jds == 1

  end function latlon_mesh_has_south_pole

  pure logical function latlon_mesh_has_north_pole(this) result(res)

    class(latlon_mesh_type), intent(in) :: this

    res = this%full_jde == global_mesh%full_nlat

  end function latlon_mesh_has_north_pole

  pure logical function latlon_mesh_is_south_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    ! FIXME: has_south_pole should be removed.
    res = j == 1

  end function latlon_mesh_is_south_pole

  pure logical function latlon_mesh_is_north_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == global_mesh%full_nlat

  end function latlon_mesh_is_north_pole

  pure logical function latlon_mesh_is_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%is_south_pole(j) .or. this%is_north_pole(j)

  end function latlon_mesh_is_pole

  pure logical function latlon_mesh_is_full_lat_next_to_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == 2 .or. j == global_mesh%full_nlat - 1

  end function latlon_mesh_is_full_lat_next_to_pole

  pure logical function latlon_mesh_is_half_lat_next_to_pole(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == 1 .or. j == global_mesh%half_nlat

  end function latlon_mesh_is_half_lat_next_to_pole

  pure logical function latlon_mesh_is_inside_with_halo_full_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= this%full_jms .and. j <= this%full_jme

  end function latlon_mesh_is_inside_with_halo_full_lat

  pure logical function latlon_mesh_is_inside_with_halo_half_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= this%half_jms .and. j <= this%half_jme

  end function latlon_mesh_is_inside_with_halo_half_lat

  pure logical function latlon_mesh_is_inside_pole_full_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= 2 .and. j <= this%full_nlat - 1

  end function latlon_mesh_is_inside_pole_full_lat

  pure logical function latlon_mesh_is_inside_pole_half_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= 1 .and. j <= this%half_nlat

  end function latlon_mesh_is_inside_pole_half_lat

  pure logical function latlon_mesh_is_outside_pole_full_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j < 1 .or. j > global_mesh%full_nlat

  end function latlon_mesh_is_outside_pole_full_lat

  pure logical function latlon_mesh_is_outside_pole_half_lat(this, j) result(res)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j < 1 .or. j > global_mesh%half_nlat

  end function latlon_mesh_is_outside_pole_half_lat

  subroutine latlon_mesh_clear(this, keep_lev)

    class(latlon_mesh_type), intent(inout) :: this
    logical, intent(in), optional :: keep_lev

    logical keep_lev_opt

    if (present(keep_lev)) then
      keep_lev_opt = keep_lev
    else
      keep_lev_opt = .false.
    end if

    if (.not. keep_lev_opt) then
      if (allocated(this%full_dlev      )) deallocate(this%full_dlev      )
      if (allocated(this%half_dlev      )) deallocate(this%half_dlev      )
      if (allocated(this%half_dlev_upper)) deallocate(this%half_dlev_upper)
      if (allocated(this%half_dlev_lower)) deallocate(this%half_dlev_lower)
      if (allocated(this%full_lev       )) deallocate(this%full_lev       )
      if (allocated(this%half_lev       )) deallocate(this%half_lev       )
    end if

    if (allocated(this%dlat            )) deallocate(this%dlat            )
    if (allocated(this%full_lon        )) deallocate(this%full_lon        )
    if (allocated(this%full_lat        )) deallocate(this%full_lat        )
    if (allocated(this%half_lon        )) deallocate(this%half_lon        )
    if (allocated(this%half_lat        )) deallocate(this%half_lat        )
    if (allocated(this%full_cos_lon    )) deallocate(this%full_cos_lon    )
    if (allocated(this%half_cos_lon    )) deallocate(this%half_cos_lon    )
    if (allocated(this%full_sin_lon    )) deallocate(this%full_sin_lon    )
    if (allocated(this%half_sin_lon    )) deallocate(this%half_sin_lon    )
    if (allocated(this%full_cos_lat    )) deallocate(this%full_cos_lat    )
    if (allocated(this%half_cos_lat    )) deallocate(this%half_cos_lat    )
    if (allocated(this%full_sin_lat    )) deallocate(this%full_sin_lat    )
    if (allocated(this%half_sin_lat    )) deallocate(this%half_sin_lat    )
    if (allocated(this%full_lon_deg    )) deallocate(this%full_lon_deg    )
    if (allocated(this%half_lon_deg    )) deallocate(this%half_lon_deg    )
    if (allocated(this%full_lat_deg    )) deallocate(this%full_lat_deg    )
    if (allocated(this%half_lat_deg    )) deallocate(this%half_lat_deg    )
    if (allocated(this%area_cell       )) deallocate(this%area_cell       )
    if (allocated(this%area_lon        )) deallocate(this%area_lon        )
    if (allocated(this%area_lon_west   )) deallocate(this%area_lon_west   )
    if (allocated(this%area_lon_east   )) deallocate(this%area_lon_east   )
    if (allocated(this%area_lon_north  )) deallocate(this%area_lon_north  )
    if (allocated(this%area_lon_south  )) deallocate(this%area_lon_south  )
    if (allocated(this%area_lat        )) deallocate(this%area_lat        )
    if (allocated(this%area_lat_west   )) deallocate(this%area_lat_west   )
    if (allocated(this%area_lat_east   )) deallocate(this%area_lat_east   )
    if (allocated(this%area_lat_north  )) deallocate(this%area_lat_north  )
    if (allocated(this%area_lat_south  )) deallocate(this%area_lat_south  )
    if (allocated(this%area_vtx        )) deallocate(this%area_vtx        )
    if (allocated(this%area_subcell    )) deallocate(this%area_subcell    )
    if (allocated(this%de_lon          )) deallocate(this%de_lon          )
    if (allocated(this%de_lat          )) deallocate(this%de_lat          )
    if (allocated(this%le_lat          )) deallocate(this%le_lat          )
    if (allocated(this%le_lon          )) deallocate(this%le_lon          )

  end subroutine latlon_mesh_clear

  subroutine latlon_mesh_final(this)

    type(latlon_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_mesh_final

end module latlon_mesh_mod
