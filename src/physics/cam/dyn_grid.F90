module dyn_grid

  use const_mod       , only: rad, radius
  use shr_kind_mod    , only: r8 => shr_kind_r8
  use pmgrid          , only: plon, plat, plev, plevp
  use cam_grid_support, only: horiz_coord_t, iMap, horiz_coord_create, &
                              cam_grid_register, cam_grid_attribute_register
  use flogger
  use string
  use const_mod       , only: p0
  use block_mod       , only: blocks, global_mesh
  use process_mod     , only: proc

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
  public get_block_levels_cnt_d
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

    type(horiz_coord_t), pointer :: lat_coord
    type(horiz_coord_t), pointer :: lon_coord
    integer(iMap), pointer :: grid_map(:,:)
    integer is, ie, js, je, i, j, ind

    plon  = global_mesh%full_nlon
    plat  = global_mesh%full_nlat
    plev  = global_mesh%full_nlev
    plevp = global_mesh%half_nlev

    lat_coord => horiz_coord_create('lat', '', plat, 'latitude',  'degree_north', 1, plat, global_mesh%full_lat_deg(1:plat))
    lon_coord => horiz_coord_create('lon', '', plon, 'longitude', 'degree_east', 1, plon, global_mesh%full_lon_deg(1:plon))

    is = blocks(1)%mesh%full_ids
    ie = blocks(1)%mesh%full_ide
    js = blocks(1)%mesh%full_jds
    je = blocks(1)%mesh%full_jde
    allocate(grid_map(4,(ie-is+1)*(je-js+1)))
    ind = 0
    do j = js, je
       do i = is, ie
          ind = ind + 1
          grid_map(1,ind) = i
          grid_map(2,ind) = j
          grid_map(3,ind) = i
          grid_map(4,ind) = j
       end do
    end do

    call cam_grid_register('latlon_grid', dyn_decomp, lat_coord, lon_coord, grid_map, unstruct=.false.)

    call cam_grid_attribute_register('latlon_grid', 'placeholder', 'what is this for?', 0)

    deallocate(grid_map)

  end subroutine dyn_grid_init

  subroutine dyn_grid_final()

  end subroutine dyn_grid_final

  subroutine dyn_grid_get_colndx(igcol, nclosest, owners, indx, jndx)

    integer, intent(in ) :: nclosest
    integer, intent(in ) :: igcol (nclosest)
    integer, intent(out) :: owners(nclosest)
    integer, intent(out) :: indx  (nclosest)
    integer, intent(out) :: jndx  (nclosest)

    stop 'dyn_grid_get_colndx not implemented'

  end subroutine dyn_grid_get_colndx

  subroutine dyn_grid_get_elem_coords( latndx, rlon, rlat, cdex )

    integer, intent(in) :: latndx ! lat  index

    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the latndx slice
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the latndx slice
    integer, optional, intent(out) :: cdex(:) ! global column index

    stop 'dyn_grid_get_elem_coords not implemented'

  end subroutine dyn_grid_get_elem_coords

  subroutine get_block_bounds_d(block_first, block_last)

    integer, intent(out) :: block_first
    integer, intent(out) :: block_last

    block_first = 1
    block_last = proc%np

  end subroutine get_block_bounds_d

  subroutine get_block_gcol_d(blockid, size, cdex)

    integer, intent(in ) :: blockid         ! Global block ID
    integer, intent(in ) :: size            ! Array size
    integer, intent(out) :: cdex(size)      ! Global column indices

    if (size /= count(proc%grid_proc_idmap == blockid)) then
      call log_error('get_block_gcol_d: size /= count(proc%grid_proc_idmap == blockid)', __FILE__, __LINE__)
    end if
    cdex = pack(proc%global_grid_id, proc%grid_proc_idmap == blockid)

  end subroutine get_block_gcol_d

  integer function get_block_gcol_cnt_d(blockid) result(res)

    integer, intent(in) :: blockid          ! Global block ID

    res = count(proc%grid_proc_idmap == blockid)

  end function get_block_gcol_cnt_d

  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

    integer, intent(in ) :: blockid         ! Global block ID
    integer, intent(in ) :: bcid            ! Column index within block
    integer, intent(in ) :: lvlsiz          ! Dimension of levels array
    integer, intent(out) :: levels(lvlsiz)  ! Levels indices for block

    integer k

    if (lvlsiz < plev + 1) then
      call log_error('get_block_levels_d: lvlsiz < plev + 1', __FILE__, __LINE__)
   else
      do k = 0, plev
         levels(k+1) = k
      end do
      do k = plev+2, lvlsiz
         levels(k) = -1
      end do
   end if

  end subroutine get_block_levels_d

  integer function get_block_levels_cnt_d(blockid, bcid) result(res)

    integer, intent(in) :: blockid          ! Global block ID
    integer, intent(in) :: bcid             ! Column index within block

    res = plevp

  end function get_block_levels_cnt_d

  integer function get_block_owner_d(blockid) result(res)

    integer, intent(in) :: blockid          ! Global block ID

    res = blockid - 1

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

    integer i, j

    j = (gcol - 1) / plon + 1
    i = gcol - (j - 1) * plon
    blockid(1) = proc%grid_proc_idmap(i,j)
    bcid(1) = proc%local_grid_id(i,j)
    do i = 2, cnt
      blockid(i) = -1
      bcid(i) = -1
    end do
    if (present(localblockid)) localblockid = -1

  end subroutine get_gcol_block_d

  integer function get_gcol_block_cnt_d(gcol) result(res)

    integer, intent(in) :: gcol             ! Global column index

    res = 1 ! One column can only be in one block.

  end function get_gcol_block_cnt_d

  subroutine get_horiz_grid_d(size, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

    integer , intent(in ) :: size                       ! Array sizes
    real(r8), intent(out), optional :: clat_d_out(size) ! Column latitudes (rad)
    real(r8), intent(out), optional :: clon_d_out(size) ! Column longitudes (rad)
    real(r8), intent(out), optional :: area_d_out(size) ! Column surface
    real(r8), intent(out), optional :: wght_d_out(size) ! Column integration
    real(r8), intent(out), optional ::  lat_d_out(size) ! Column latitudes (deg)
    real(r8), intent(out), optional ::  lon_d_out(size) ! Column longitudes (deg)

    integer i, j, icol

    if (present(clat_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          clat_d_out(icol) = global_mesh%full_lat(j)
          icol = icol + 1
        end do
      end do
    end if
    if (present(clon_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          clon_d_out(icol) = global_mesh%full_lon(i)
          icol = icol + 1
        end do
      end do
    end if
    if (present(area_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          area_d_out(icol) = global_mesh%area_cell(j) / radius**2
          icol = icol + 1
        end do
      end do
    end if
    if (present(wght_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          wght_d_out(icol) = global_mesh%area_cell(j) / radius**2
          icol = icol + 1
        end do
      end do
    end if
    if (present( lat_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          lat_d_out(icol) = global_mesh%full_lat_deg(j)
          icol = icol + 1
        end do
      end do
    end if
    if (present( lon_d_out)) then
      icol = 1
      do j = global_mesh%full_jds, global_mesh%full_jde
        do i = global_mesh%full_ids, global_mesh%full_ide
          lon_d_out(icol) = global_mesh%full_lon_deg(i)
          icol = icol + 1
        end do
      end do
    end if

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
    grid_attribute_names(1) = 'placeholder'

  end subroutine physgrid_copy_attributes_d

  function get_dyn_grid_parm_real1d(name) result(rval)

    character(*), intent(in) :: name
    real(r8), pointer :: rval(:)

    stop 'get_dyn_grid_parm_real1d not implemented'
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