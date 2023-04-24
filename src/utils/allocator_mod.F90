module allocator_mod

  use flogger
  use mesh_mod

  implicit none

  private

  public allocate_array

  interface allocate_array
    module procedure allocate_array_1d_r4
    module procedure allocate_array_1d_r4_extra_dim
    module procedure allocate_array_2d_r4
    module procedure allocate_array_2d_r4_extra_dim
    module procedure allocate_array_3d_r4
    module procedure allocate_array_3d_r4_extra_dim
    module procedure allocate_array_1d_r8
    module procedure allocate_array_1d_r8_extra_dim
    module procedure allocate_array_2d_r8
    module procedure allocate_array_2d_r8_extra_dim
    module procedure allocate_array_3d_r8
    module procedure allocate_array_3d_r8_extra_dim
    module procedure allocate_array_3d_i4
    module procedure allocate_array_pointer_3d_r4
    module procedure allocate_array_pointer_3d_r8
    module procedure allocate_array_1d_r16
    module procedure allocate_array_2d_r16
    module procedure allocate_array_3d_r16
  end interface allocate_array

contains

  subroutine allocate_array_1d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_ims:mesh%full_ime))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_ims:mesh%half_ime))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_jms:mesh%full_jme))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_jms:mesh%half_jme))
    else
      call log_error('allocate_array_1d_r4: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r4

  subroutine allocate_array_1d_r4_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_ims:mesh%full_ime, extra_dim))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_ims:mesh%half_ime, extra_dim))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_jms:mesh%half_jme, extra_dim))
    else
      call log_error('allocate_array_1d_r4_extra_dim: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r4_extra_dim

  subroutine allocate_array_2d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r4

  subroutine allocate_array_2d_r4_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme, extra_dim))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme, extra_dim))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r4_extra_dim

  subroutine allocate_array_3d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r4

  subroutine allocate_array_3d_r4_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), allocatable :: array(:,:,:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r4_extra_dim

  subroutine allocate_array_1d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_ims:mesh%full_ime))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_ims:mesh%half_ime))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_jms:mesh%full_jme))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_jms:mesh%half_jme))
    else
      call log_error('allocate_array_1d_r8: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r8

  subroutine allocate_array_1d_r8_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_ims:mesh%full_ime, extra_dim))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_ims:mesh%half_ime, extra_dim))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_jms:mesh%half_jme, extra_dim))
    else
      call log_error('allocate_array_1d_r8_extra_dim: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r8_extra_dim

  subroutine allocate_array_2d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r8

  subroutine allocate_array_2d_r8_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme, extra_dim))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme, extra_dim))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme, extra_dim))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r8_extra_dim

  subroutine allocate_array_3d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r8

  subroutine allocate_array_3d_r8_extra_dim(mesh, array, extra_dim, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), allocatable :: array(:,:,:,:)
    integer        , intent(in )              :: extra_dim
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme, extra_dim))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme, extra_dim))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r8_extra_dim

  subroutine allocate_array_3d_i4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    integer(4)     , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0

  end subroutine allocate_array_3d_i4

  subroutine allocate_array_pointer_3d_r4(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(4)        , intent(out), pointer     :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0.0d0

  end subroutine allocate_array_pointer_3d_r4

  subroutine allocate_array_pointer_3d_r8(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(8)        , intent(out), pointer     :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0.0d0

  end subroutine allocate_array_pointer_3d_r8

  subroutine allocate_array_1d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (allocated(array)) deallocate(array)

    if (present(full_lon)) then
      if (full_lon) allocate(array(mesh%full_ims:mesh%full_ime))
    else if (present(half_lon)) then
      if (half_lon) allocate(array(mesh%half_ims:mesh%half_ime))
    else if (present(full_lat)) then
      if (full_lat) allocate(array(mesh%full_jms:mesh%full_jme))
    else if (present(half_lat)) then
      if (half_lat) allocate(array(mesh%half_jms:mesh%half_jme))
    else
      call log_error('allocate_array_1d_r16: Missing full_lon, half_lon, full_lat or half_lat argument!')
    end if

    array = 0.0d0

  end subroutine allocate_array_1d_r16

  subroutine allocate_array_2d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat

    if (present(full_lon) .and. present(full_lat)) then
      if (full_lon .and. full_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    else if (present(full_lon) .and. present(half_lat)) then
      if (full_lon .and. half_lat) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme))
    else if (present(half_lon) .and. present(full_lat)) then
      if (half_lon .and. full_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme))
    else if (present(half_lon) .and. present(half_lat)) then
      if (half_lon .and. half_lat) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
    end if

    array = 0.0d0

  end subroutine allocate_array_2d_r16

  subroutine allocate_array_3d_r16(mesh, array, full_lon, half_lon, full_lat, half_lat, full_lev, half_lev)

    type(mesh_type), intent(in )              :: mesh
    real(16)        , intent(out), allocatable :: array(:,:,:)
    logical        , intent(in ), optional    :: full_lon
    logical        , intent(in ), optional    :: half_lon
    logical        , intent(in ), optional    :: full_lat
    logical        , intent(in ), optional    :: half_lat
    logical        , intent(in ), optional    :: full_lev
    logical        , intent(in ), optional    :: half_lev

    if (present(full_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (full_lon .and. full_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (full_lon .and. half_lat .and. full_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(full_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (full_lon .and. full_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(full_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (full_lon .and. half_lat .and. half_lev) allocate(array(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(full_lev)) then
      if (half_lon .and. full_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(full_lev)) then
      if (half_lon .and. half_lat .and. full_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
    else if (present(half_lon) .and. present(full_lat) .and. present(half_lev)) then
      if (half_lon .and. full_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
    else if (present(half_lon) .and. present(half_lat) .and. present(half_lev)) then
      if (half_lon .and. half_lat .and. half_lev) allocate(array(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme,mesh%half_kms:mesh%half_kme))
    else
      allocate(array(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
    end if

    array = 0.0d0

  end subroutine allocate_array_3d_r16

end module allocator_mod
