module dyn_grid

  use const_mod, only: rad
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plon, plat, plev, plevp
  use hycoef, only: hycoef_init, hycoef_final
  use const_mod, only: p0
  use block_mod
  use vert_coord_mod, only: hyai, hybi

  implicit none

  private

  public dyn_grid_init
  public dyn_grid_final
  public dyn_grid_get_colndx
  public dyn_grid_get_elem_coords
  public get_block_bounds_d
  public get_block_gcol_d
  public get_block_gcol_cnt_d
  public get_block_levels_d
  public get_block_lvl_cnt_d
  public get_block_owner_d
  public get_dyn_grid_parm
  public get_gcol_block_d
  public get_gcol_block_cnt_d
  public get_horiz_grid_d
  public get_horiz_grid_dim_d
  public physgrid_copy_attributes_d
  public get_dyn_grid_parm_real1d

  integer, parameter, public :: dyn_decomp    = 101
  integer, parameter, public :: ptimelevels   = 2     ! number of time levels in the dycore

contains

  subroutine dyn_grid_init()

    plon  = blocks(1)%mesh%full_nlon
    plat  = blocks(1)%mesh%full_nlat
    plev  = blocks(1)%mesh%full_nlev
    plevp = blocks(1)%mesh%half_nlev

    call hycoef_init(hyai, hybi, p0, p0)

  end subroutine dyn_grid_init

  subroutine dyn_grid_final()

    call hycoef_final()

  end subroutine dyn_grid_final

  subroutine dyn_grid_get_colndx(igcol, nclosest, owners, indx, jndx)
  
    integer, intent(in ) :: nclosest
    integer, intent(in ) :: igcol (nclosest)
    integer, intent(out) :: owners(nclosest)
    integer, intent(out) :: indx  (nclosest)
    integer, intent(out) :: jndx  (nclosest)

  end subroutine dyn_grid_get_colndx

  subroutine dyn_grid_get_elem_coords( latndx, rlon, rlat, cdex )
  
    integer, intent(in) :: latndx ! lat  index
  
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the latndx slice
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the latndx slice
    integer, optional, intent(out) :: cdex(:) ! global column index

  end subroutine dyn_grid_get_elem_coords

  subroutine get_block_bounds_d(block_first, block_last)

    integer, intent(out) :: block_first
    integer, intent(out) :: block_last

  end subroutine get_block_bounds_d

  subroutine get_block_gcol_d(blockid, size, cdex)

    integer, intent(in ) :: blockid         ! Global block ID
    integer, intent(in ) :: size            ! Array size
    integer, intent(out) :: cdex(size)      ! Global column indices

  end subroutine get_block_gcol_d

  integer function get_block_gcol_cnt_d(blockid) result(res)

    integer, intent(in) :: blockid          ! Global block ID

  end function get_block_gcol_cnt_d

  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

    integer, intent(in ) :: blockid         ! Global block ID
    integer, intent(in ) :: bcid            ! Column index within block
    integer, intent(in ) :: lvlsiz          ! Dimension of levels array
    integer, intent(out) :: levels(lvlsiz)  ! Levels indices for block

  end subroutine get_block_levels_d

  integer function get_block_lvl_cnt_d(blockid, bcid) result(res)

    integer, intent(in) :: blockid          ! Global block ID
    integer, intent(in) :: bcid             ! Column index within block

    res = plevp

  end function get_block_lvl_cnt_d

  integer function get_block_owner_d(blockid) result(res)

    integer, intent(in) :: blockid          ! Global block ID

  end function get_block_owner_d

  integer function get_dyn_grid_parm(name) result(res)

    character(*), intent(in) :: name

    select case (name)
    case ('plon')
      res = plon
    case ('plat')
      res = plat
    case ('plev')
      res = plev
    case ('plevp')
      res = plevp
    case ('beglonxy')
      res = blocks(1)%mesh%full_ids
    case ('endlonxy')
      res = blocks(1)%mesh%full_ide
    case ('beglatxy')
      res = blocks(1)%mesh%full_jds
    case ('endlatxy')
      res = blocks(1)%mesh%full_jde
    case default
      res = -1
    end select

  end function get_dyn_grid_parm

  subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)

    ! Purpose: Return global block index and local column index for global column index.

    integer, intent(in ) :: gcol            ! Global column index
    integer, intent(in ) :: cnt             ! Size of blockid and bcid arrays
    integer, intent(out) :: blockid(cnt)    ! Block index
    integer, intent(out) :: bcid   (cnt)    ! Column index within block
    integer, intent(out), optional :: localblockid(cnt)

  end subroutine get_gcol_block_d

  integer function get_gcol_block_cnt_d(gcol) result(res)

    integer, intent(in) :: gcol             ! Global column index

    res = 1

  end function get_gcol_block_cnt_d

  subroutine get_horiz_grid_d(size, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

    integer , intent(in ) :: size                       ! Array sizes
    real(r8), intent(out), optional :: clat_d_out(size) ! Column latitudes
    real(r8), intent(out), optional :: clon_d_out(size) ! Column longitudes
    real(r8), intent(out), optional :: area_d_out(size) ! Column surface
    real(r8), intent(out), optional :: wght_d_out(size) ! Column integration
    real(r8), intent(out), optional :: lat_d_out(size)  ! Column deg latitudes
    real(r8), intent(out), optional :: lon_d_out(size)  ! Column deg longitudes

  end subroutine get_horiz_grid_d

  subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)

    integer, intent(out) :: hdim1_d
    integer, intent(out) :: hdim2_d

    hdim1_d = plon
    hdim2_d = plat

  end subroutine get_horiz_grid_dim_d

  subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)

    character(*), intent(out) :: gridname
    character(*), intent(inout), pointer :: grid_attribute_names(:)

    gridname = 'latlon_grid'
    allocate(grid_attribute_names(1))
    grid_attribute_names(1) = 'gw'

  end subroutine physgrid_copy_attributes_d

  function get_dyn_grid_parm_real1d(name) result(rval)

    character(*), intent(in) :: name
    real(r8), pointer :: rval(:)

    ! if (name == 'clat') then
    !   rval => clat
    ! else if (name == 'latdeg') then
    !   rval => latdeg
    ! else if (name == 'w') then
    !   rval => w
    ! else
    !   nullify(rval)
    ! end if

  end function get_dyn_grid_parm_real1d

end module dyn_grid