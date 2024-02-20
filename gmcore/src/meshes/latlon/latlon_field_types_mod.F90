! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This is module provides field types which encapsulates meta data with array.
!   By doing so, users can send the field objects into their functions without
!   the need to specify the array index ranges.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module latlon_field_types_mod

  use fiona
  use flogger
  use const_mod
  use latlon_mesh_mod
  use latlon_halo_mod

  implicit none

  private

  public latlon_mesh_type
  public latlon_halo_type
  public latlon_field2d_type
  public latlon_field3d_type
  public latlon_field4d_type

  integer, public, parameter :: field_name_len      = 32
  integer, public, parameter :: field_long_name_len = 128
  integer, public, parameter :: field_units_len     = 32
  integer, public, parameter :: field_loc_len       = 10

  type latlon_field_meta_type
    character(field_name_len     ) :: name      = 'N/A'
    character(field_long_name_len) :: long_name = 'N/A'
    character(field_units_len    ) :: units     = 'N/A'
    character(field_loc_len      ) :: loc       = 'N/A'
    logical :: full_lon         = .true.
    logical :: full_lat         = .true.
    logical :: full_lev         = .true.
    logical :: initialized      = .false.
    logical :: linked           = .false.
    logical :: restart          = .false.
    logical :: initial          = .false.
    logical :: halo_cross_pole  = .false.
    type(latlon_mesh_type), pointer :: mesh     => null()
    type(latlon_halo_type), pointer :: halo(:)  => null()
  end type

  type, extends(latlon_field_meta_type) :: latlon_field2d_type
    real(r8), contiguous, pointer :: d(:,:) => null()
  contains
    procedure :: init  => latlon_field2d_init
    procedure :: clear => latlon_field2d_clear
    procedure, private :: latlon_field2d_link_2d
    procedure, private :: latlon_field2d_link_3d
    generic :: link => latlon_field2d_link_2d, latlon_field2d_link_3d
    final latlon_field2d_final
  end type latlon_field2d_type

  type, extends(latlon_field_meta_type) :: latlon_field3d_type
    real(r8), contiguous, pointer :: d(:,:,:) => null()
  contains
    procedure :: init  => latlon_field3d_init
    procedure :: clear => latlon_field3d_clear
    procedure, private :: latlon_field3d_link_3d
    procedure, private :: latlon_field3d_link_4d
    generic :: link => latlon_field3d_link_3d, latlon_field3d_link_4d
    final latlon_field3d_final
  end type latlon_field3d_type

  type, extends(latlon_field_meta_type) :: latlon_field4d_type
    real(r8), contiguous, pointer :: d(:,:,:,:) => null()
  contains
    procedure :: init  => latlon_field4d_init
    procedure :: clear => latlon_field4d_clear
    procedure :: link  => latlon_field4d_link
    final latlon_field4d_final
  end type latlon_field4d_type

contains

  subroutine latlon_field2d_init(this, name, long_name, units, loc, mesh, halo, halo_cross_pole)

    class(latlon_field2d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole

    call this%clear()

    this%name      = name
    this%long_name = long_name
    this%units     = units
    this%loc       = loc
    this%mesh      => mesh
    this%halo      => halo
    if (present(halo_cross_pole)) this%halo_cross_pole = halo_cross_pole

    select case (loc)
    case ('cell')
      this%full_lon = .true. ; this%full_lat = .true.
      allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    case ('lon')
      this%full_lon = .false.; this%full_lat = .true.
      allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme))
    case ('lat')
      this%full_lon = .true. ; this%full_lat = .false.
      allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme))
    case ('vtx')
      this%full_lon = .false.; this%full_lat = .false.
      allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme))
    end select

    this%d = 0
    this%initialized = .true.

  end subroutine latlon_field2d_init

  subroutine latlon_field2d_clear(this)

    class(latlon_field2d_type), intent(inout) :: this

    if (this%initialized .and. .not. this%linked .and. associated(this%d)) then
      deallocate(this%d)
      this%d => null()
    end if
    this%name            = 'N/A'
    this%long_name       = 'N/A'
    this%units           = 'N/A'
    this%loc             = 'N/A'
    this%mesh            => null()
    this%halo            => null()
    this%halo_cross_pole = .false.
    this%full_lon        = .true.
    this%full_lat        = .true.
    this%initialized     = .false.
    this%linked          = .false.
    this%restart         = .false.
    this%initial         = .false.

  end subroutine latlon_field2d_clear

  subroutine latlon_field2d_link_2d(this, other)

    class(latlon_field2d_type), intent(inout) :: this
    type(latlon_field2d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field2d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      this%name            = other%name
      this%long_name       = other%long_name
      this%units           = other%units
      this%loc             = other%loc
      this%mesh            => other%mesh
      this%halo            => other%halo
      this%halo_cross_pole = other%halo_cross_pole
      this%full_lon        = other%full_lon
      this%full_lat        = other%full_lat
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d => other%d
    this%linked = .true.
    this%initialized = .true.

  end subroutine latlon_field2d_link_2d

  subroutine latlon_field2d_link_3d(this, other, i3)

    class(latlon_field2d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other
    integer, intent(in) :: i3

    real(r8), pointer, contiguous :: tmp(:,:)
    integer is, ie, js, je

    if (this%initialized .and. .not. (this%full_lon .eqv. other%full_lon .and. this%full_lat .eqv. other%full_lat)) then
      call log_error('latlon_field2d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    select case (this%loc)
    case ('cell')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
    case ('lon')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
    case ('lat')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
    case ('vtx')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
    end select
    ! Use a temporary array pointer to fix compile error.
    tmp => other%d(:,:,i3)
    this%d(is:ie,js:je) => tmp
    this%linked = .true.
    this%initialized = .true.

  end subroutine latlon_field2d_link_3d

  subroutine latlon_field2d_final(this)

    type(latlon_field2d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field2d_final

  subroutine latlon_field3d_init(this, name, long_name, units, loc, mesh, halo, halo_cross_pole, ptr_to)

    class(latlon_field3d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    type(latlon_field3d_type), intent(in), optional, target :: ptr_to

    call this%clear()

    this%name      = name
    this%long_name = long_name
    this%units     = units
    this%loc       = loc
    this%mesh      => mesh
    this%halo      => halo
    if (present(halo_cross_pole)) this%halo_cross_pole = halo_cross_pole

    if (present(ptr_to)) then
      this%d => ptr_to%d
      this%linked = .true.
    else
      select case (loc)
      case ('cell')
        this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .true.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      case ('lon')
        this%full_lon = .false.; this%full_lat = .true. ; this%full_lev = .true.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      case ('lat')
        this%full_lon = .true. ; this%full_lat = .false.; this%full_lev = .true.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
      case ('vtx')
        this%full_lon = .false.; this%full_lat = .false.; this%full_lev = .true.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
      case ('lev')
        this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .false.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
      case ('lev_lon')
        this%full_lon = .false.; this%full_lat = .true. ; this%full_lev = .false.
        allocate(this%d(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
      case ('lev_lat')
        this%full_lon = .true. ; this%full_lat = .false.; this%full_lev = .false.
        allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
      end select
    end if

    this%d = 0
    this%initialized = .true.

  end subroutine latlon_field3d_init

  subroutine latlon_field3d_clear(this)

    class(latlon_field3d_type), intent(inout) :: this

    if (this%initialized .and. .not. this%linked .and. associated(this%d)) then
      deallocate(this%d)
      this%d => null()
    end if
    this%name            = 'N/A'
    this%long_name       = 'N/A'
    this%units           = 'N/A'
    this%loc             = 'N/A'
    this%mesh            => null()
    this%halo            => null()
    this%halo_cross_pole = .false.
    this%full_lon        = .true.
    this%full_lat        = .true.
    this%full_lev        = .true.
    this%initialized     = .false.
    this%linked          = .false.
    this%restart         = .false.
    this%initial         = .false.

  end subroutine latlon_field3d_clear

  subroutine latlon_field3d_link_3d(this, other)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field3d_link_3d: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      this%name            = other%name
      this%long_name       = other%long_name
      this%units           = other%units
      this%loc             = other%loc
      this%mesh            => other%mesh
      this%halo            => other%halo
      this%halo_cross_pole = other%halo_cross_pole
      this%full_lon        = other%full_lon
      this%full_lat        = other%full_lat
      this%full_lev        = other%full_lev
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d => other%d
    this%linked = .true.
    this%initialized = .true.

  end subroutine latlon_field3d_link_3d

  subroutine latlon_field3d_link_4d(this, other, i4)

    class(latlon_field3d_type), intent(inout) :: this
    type(latlon_field4d_type), intent(in) :: other
    integer, intent(in) :: i4

    real(r8), pointer, contiguous :: tmp(:,:,:)
    integer is, ie, js, je, ks, ke

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field3d_link_4d: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      this%name            = other%name
      this%long_name       = other%long_name
      this%units           = other%units
      this%loc             = other%loc
      this%mesh            => other%mesh
      this%halo            => other%halo
      this%halo_cross_pole = other%halo_cross_pole
      this%full_lon        = other%full_lon
      this%full_lat        = other%full_lat
      this%full_lev        = other%full_lev
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    select case (this%loc)
    case ('cell')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lon')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lat')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    case ('lev')
      is = this%mesh%full_ims; ie = this%mesh%full_ime
      js = this%mesh%full_jms; je = this%mesh%full_jme
      ks = this%mesh%half_kms; ke = this%mesh%half_kme
    case ('vtx')
      is = this%mesh%half_ims; ie = this%mesh%half_ime
      js = this%mesh%half_jms; je = this%mesh%half_jme
      ks = this%mesh%full_kms; ke = this%mesh%full_kme
    end select
    ! Use a temporary array pointer to fix compile error.
    tmp => other%d(:,:,:,i4)
    this%d(is:ie,js:je,ks:ke) => tmp
    this%linked = .true.
    this%initialized = .true.

  end subroutine latlon_field3d_link_4d

  subroutine latlon_field3d_final(this)

    type(latlon_field3d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field3d_final

  subroutine latlon_field4d_init(this, name, long_name, units, loc, mesh, halo, halo_cross_pole, n4)

    class(latlon_field4d_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in), target :: halo(:)
    logical, intent(in), optional :: halo_cross_pole
    integer, intent(in) :: n4

    call this%clear()

    this%name      = name
    this%long_name = long_name
    this%units     = units
    this%loc       = loc
    this%mesh      => mesh
    this%halo      => halo
    if (present(halo_cross_pole)) this%halo_cross_pole = halo_cross_pole

    select case (loc)
    case ('cell')
      this%full_lon = .true. ; this%full_lat = .true. ; this%full_lev = .true.
      allocate(this%d(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme,n4))
    end select

    this%d = 0
    this%linked = .false.
    this%initialized = .true.

  end subroutine latlon_field4d_init

  subroutine latlon_field4d_clear(this)

    class(latlon_field4d_type), intent(inout) :: this

    if (this%initialized .and. .not. this%linked .and. associated(this%d)) then
      deallocate(this%d)
      this%d => null()
    end if
    this%name            = 'N/A'
    this%long_name       = 'N/A'
    this%units           = 'N/A'
    this%loc             = 'N/A'
    this%mesh            => null()
    this%halo            => null()
    this%halo_cross_pole = .false.
    this%initialized     = .false.
    this%linked          = .false.
    this%restart         = .false.
    this%initial         = .false.

  end subroutine latlon_field4d_clear

  subroutine latlon_field4d_link(this, other)

    class(latlon_field4d_type), intent(inout) :: this
    type(latlon_field4d_type), intent(in) :: other

    if (this%initialized .and. this%loc /= other%loc) then
      call log_error('latlon_field4d_link: cannot link fields with different loc!', __FILE__, __LINE__)
    else
      this%name            = other%name
      this%long_name       = other%long_name
      this%units           = other%units
      this%loc             = other%loc
      this%mesh            => other%mesh
      this%halo            => other%halo
      this%halo_cross_pole = other%halo_cross_pole
      this%full_lon        = other%full_lon
      this%full_lat        = other%full_lat
      this%full_lev        = other%full_lev
    end if
    if (this%initialized .and. .not. this%linked .and. associated(this%d)) deallocate(this%d)
    this%d => other%d
    this%linked = .true.
    this%initialized = .true.

  end subroutine latlon_field4d_link

  subroutine latlon_field4d_final(this)

    type(latlon_field4d_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_field4d_final

end module latlon_field_types_mod
