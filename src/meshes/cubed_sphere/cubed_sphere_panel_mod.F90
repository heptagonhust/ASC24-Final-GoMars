! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module cubed_sphere_panel_mod

  use const_mod
  use math_mod
  use sphere_geometry_mod

  implicit none

  private

  public cubed_sphere_panel_type

  type ngb_type
    type(cubed_sphere_panel_type), pointer :: ptr => null()
  end type ngb_type

  type grid_type
    real(r8), allocatable, dimension(  :,:,:) :: p        ! Grid point vector (i.e. Cartesian coordinate)
    real(r8), allocatable, dimension(      :) :: x1       ! x on cube tile for panel 1
    real(r8), allocatable, dimension(      :) :: y1       ! y on cube tile for panel 1
    real(r8), allocatable, dimension(    :,:) :: lx       ! Edge length along x axis
    real(r8), allocatable, dimension(    :,:) :: ly       ! Edge length along y axis
    real(r8), allocatable, dimension(    :,:) :: lon      ! Longitude (radian)
    real(r8), allocatable, dimension(    :,:) :: lon_deg  ! Longitude (degree)
    real(r8), allocatable, dimension(    :,:) :: lat      ! Latitude (radian)
    real(r8), allocatable, dimension(    :,:) :: lat_deg  ! Latitude (degree)
    real(r8), allocatable, dimension(    :,:) :: area     ! Area of cell
    real(r8), allocatable, dimension(  :,:,:) :: subarea  ! Subarea of cell
    real(r8), allocatable, dimension(    :,:) :: sina     ! Sine of angle between local unit base vectors
    real(r8), allocatable, dimension(    :,:) :: cosa     ! Cosine of angle between local unit base vectors
    real(r8), allocatable, dimension(:,:,:,:) :: A        ! Transform matrix from cube to sphere
    real(r8), allocatable, dimension(:,:,:,:) :: iA       ! Transform matrix from sphere to cube
  contains
    procedure :: init         => grid_init
    procedure :: clear        => grid_clear
    procedure :: set_lon_lat  => grid_set_lon_lat
    final :: grid_final
  end type grid_type

  type cubed_sphere_panel_type
    logical :: initialized  = .false.
    logical :: active       = .false.
    integer :: id           = 0
    integer :: proj_type    = 1
    integer :: nx  = 0, ny  = 0, nz  = 0
    integer :: hw  = 0
    integer :: full_ids = 0, full_ide = 0, half_ids = 0, half_ide = 0
    integer :: full_jds = 0, full_jde = 0, half_jds = 0, half_jde = 0
    integer :: full_kds = 0, full_kde = 0, half_kds = 0, half_kde = 0
    integer :: full_ndx = 0, full_ndy = 0, half_ndx = 0, half_ndy = 0
    integer :: full_ims = 0, full_ime = 0, half_ims = 0, half_ime = 0
    integer :: full_jms = 0, full_jme = 0, half_jms = 0, half_jme = 0
    integer :: full_kms = 0, full_kme = 0, half_kms = 0, half_kme = 0
    integer :: full_nmx = 0, full_nmy = 0, half_nmx = 0, half_nmy = 0
    type(grid_type) cell, vtx, xedge, yedge
    type(ngb_type) ngb(4)
  contains
    procedure :: init    => cubed_sphere_panel_init
    procedure :: clear   => cubed_sphere_panel_clear
    procedure :: connect => cubed_sphere_panel_connect
    final :: cubed_sphere_panel_final
  end type cubed_sphere_panel_type

contains

  recursive subroutine cubed_sphere_panel_init(this, proj_type, nc, nz, hw, ids, jds, ndx, ndy, active)

    class(cubed_sphere_panel_type), intent(inout) :: this
    integer, intent(in) :: proj_type      ! Projection type: 1 - equiangular; 2 - equi-edge
    integer, intent(in) :: nc             ! Panel grid size along x axis and y axis
    integer, intent(in) :: nz             ! Level size
    integer, intent(in) :: hw             ! Halo width
    integer, intent(in), optional :: ids  ! Domain start index along x axis
    integer, intent(in), optional :: jds  ! Domain start index along y axis
    integer, intent(in), optional :: ndx  ! Domain grid size along x axis
    integer, intent(in), optional :: ndy  ! Domain grid size along y axis
    logical, intent(in), optional :: active

    real(r8) alpha, Rref, dela, e1(3), e2(3), elon(3), elat(3), G(2,2), iG(2,2), A(2,2)
    real(r8) sinlon, coslon, sinlat, coslat
    integer i, j, ngb_ids, ngb_jds, ngb_ndx, ngb_ndy
    logical even_panel

    call this%clear()

    this%proj_type = proj_type
    this%nx        = nc
    this%ny        = nc
    this%nz        = nz
    this%hw        = hw

    if (present(active)) this%active = active

    if (present(ids) .and. present(ndx)) then
      this%full_ids = ids; this%full_ide = ids + ndx - 1; this%full_ndx = ndx
    else
      this%full_ids = 1  ; this%full_ide = nc           ; this%full_ndx = nc
    end if
    this%half_ids = this%full_ids; this%half_ide = this%full_ide + 1; this%half_ndx = this%full_ndx + 1
    if (present(jds) .and. present(ndy)) then
      this%full_jds = jds; this%full_jde = jds + ndy - 1; this%full_ndy = ndy
    else
      this%full_jds = 1  ; this%full_jde = nc           ; this%full_ndy = nc
    end if
    this%half_jds = this%full_jds; this%half_jde = this%full_jde + 1; this%half_ndy = this%full_ndy + 1
    this%full_kds = 1; this%full_kde = nz
    this%half_kds = 1; this%half_kde = nz + 1

    this%full_ims = this%full_ids - hw; this%full_ime = this%full_ide + hw
    this%full_jms = this%full_jds - hw; this%full_jme = this%full_jde + hw
    this%full_kms = this%full_kds - hw; this%full_kme = this%full_kde + hw
    this%half_ims = this%half_ids - hw; this%half_ime = this%half_ide + hw
    this%half_jms = this%half_jds - hw; this%half_jme = this%half_jde + hw
    this%half_kms = this%half_kds - hw; this%half_kme = this%half_kde + hw

    this%full_nmx = this%full_ime - this%full_ims + 1
    this%full_nmy = this%full_jme - this%full_jms + 1
    this%half_nmx = this%half_ime - this%half_ims + 1
    this%half_nmy = this%half_jme - this%half_jms + 1

    ! Initialize grids on different locations.
    ! NOTE: Enlarge one level of halo for mesh quantity calculations.
    call this%cell %init(this%full_ims-1, this%full_ime+1, this%full_jms-1, this%full_jme+1)
    call this%vtx  %init(this%half_ims-1, this%half_ime+1, this%half_jms-1, this%half_jme+1)
    call this%xedge%init(this%half_ims-1, this%half_ime+1, this%full_jms-1, this%full_jme+1)
    call this%yedge%init(this%full_ims-1, this%full_ime+1, this%half_jms-1, this%half_jme+1)

    ! --------------------------------------------------------------------------
    ! Setup grid coordinates.
    select case (proj_type)
    case (1) ! Equiangular grid
      alpha = pi / 4.0_r8
      Rref  = 1.0_r8
    case (2) ! Equi-edge grid
      alpha = asin(1.0_r8 / sqrt(3.0_r8))
      Rref  = sqrt(2.0_r8)
    case default
    end select
    dela = 2 * alpha / this%nx
    ! Cell
    do i = this%full_ims - 1, this%full_ime + 1
      this%cell%x1(i) = Rref * tan(dela * (i - 0.5_r8) - alpha)
    end do
    do j = this%full_jms - 1, this%full_jme + 1
      this%cell%y1(j) = Rref * tan(dela * (j - 0.5_r8) - alpha)
    end do
    call this%cell%set_lon_lat(this%id)
    ! Vertex
    do i = this%half_ims - 1, this%half_ime + 1
      this%vtx%x1(i) = Rref * tan(dela * (i - 1.0_r8) - alpha)
    end do
    do j = this%half_jms - 1, this%half_jme + 1
      this%vtx%y1(j) = Rref * tan(dela * (j - 1.0_r8) - alpha)
    end do
    call this%vtx%set_lon_lat(this%id)
    ! X edge
    do i = this%half_ims - 1, this%half_ime + 1
      this%xedge%x1(i) = Rref * tan(dela * (i - 1.0_r8) - alpha)
    end do
    do j = this%full_jms - 1, this%full_jme + 1
      this%xedge%y1(j) = Rref * tan(dela * (j - 0.5_r8) - alpha)
    end do
    call this%xedge%set_lon_lat(this%id)
    ! Y edge
    do i = this%full_ims - 1, this%full_ime + 1
      this%yedge%x1(i) = Rref * tan(dela * (i - 0.5_r8) - alpha)
    end do
    do j = this%half_jms - 1, this%half_jme + 1
      this%yedge%y1(j) = Rref * tan(dela * (j - 1.0_r8) - alpha)
    end do
    call this%yedge%set_lon_lat(this%id)

    ! --------------------------------------------------------------------------
    ! Setup edge lengths.
    ! X edge (ly)
    do j = this%full_jms - 1, this%full_jme + 1
      do i = this%half_ims - 1, this%half_ime + 1
        this%xedge%ly(i,j) = great_circle(radius  , &
          this%vtx%lon(i,j  ), this%vtx%lat(i,j  ), &
          this%vtx%lon(i,j+1), this%vtx%lat(i,j+1))
      end do
    end do
    ! Y edge (lx)
    do j = this%half_jms - 1, this%half_jme + 1
      do i = this%full_ims - 1, this%full_ime + 1
        this%yedge%lx(i,j) = great_circle(radius  , &
          this%vtx%lon(i  ,j), this%vtx%lat(i  ,j), &
          this%vtx%lon(i+1,j), this%vtx%lat(i+1,j))
      end do
    end do
    ! Average edge lengthes to other places.
    ! Cell (lx and ly)
    do j = this%full_jms - 1, this%full_jme + 1
      do i = this%full_ims - 1, this%full_ime + 1
        this%cell%lx(i,j) = (this%yedge%lx(i,j) + this%yedge%lx(i,j+1)) * 0.5_r8
        this%cell%ly(i,j) = (this%xedge%ly(i,j) + this%xedge%ly(i+1,j)) * 0.5_r8
      end do
    end do
    ! X edge (lx)
    do j = this%full_jms, this%full_jme
      do i = this%half_ims, this%half_ime
        this%xedge%lx(i,j) = (this%cell%lx(i-1,j) + this%cell%lx(i,j)) * 0.5_r8
      end do
    end do
    ! Y edge (ly)
    do j = this%half_jms, this%half_jme
      do i = this%full_ims, this%full_ime
        this%yedge%ly(i,j) = (this%cell%ly(i,j-1) + this%cell%ly(i,j)) * 0.5_r8
      end do
    end do

    ! --------------------------------------------------------------------------
    ! Setup sine and cosine of angles between local unit base vectors, and transform matrices.
    ! Cell
    do j = this%full_jms, this%full_jme
      do i = this%full_ims, this%full_ime
        associate (p0 => this%cell %p(:,i  ,j  ), &
                   p1 => this%xedge%p(:,i  ,j  ), &
                   p2 => this%xedge%p(:,i+1,j  ), &
                   p3 => this%yedge%p(:,i  ,j  ), &
                   p4 => this%yedge%p(:,i  ,j+1))
        e1 = cross_product(cross_product(p1, p2), p0); e1 = e1 / norm2(e1)
        e2 = cross_product(cross_product(p3, p4), p0); e2 = e2 / norm2(e2)
        this%cell%sina(i,j) = dot_product(p0, cross_product(e1, e2))
        this%cell%cosa(i,j) = dot_product(e1, e2)
        sinlon = sin(this%cell%lon(i,j)); coslon = cos(this%cell%lon(i,j))
        sinlat = sin(this%cell%lat(i,j)); coslat = cos(this%cell%lat(i,j))
        G  = reshape([1.0_r8,  this%cell%cosa(i,j),  this%cell%cosa(i,j), 1.0_r8], [2,2])
        iG = reshape([1.0_r8, -this%cell%cosa(i,j), -this%cell%cosa(i,j), 1.0_r8], [2,2]) / this%cell%sina(i,j)**2
        elon = [-sinlon, coslon, 0.0_r8]
        elat = [-sinlat*coslon, -sinlat*sinlon, coslat]
        A = reshape([dot_product(e1,elon),dot_product(e1,elat),dot_product(e2,elon),dot_product(e2,elat)], [2,2])
        this%cell%A(:,:,i,j) = matmul(A, iG)
        A = reshape([this%cell%A(1,1,i,j),-this%cell%A(2,1,i,j),-this%cell%A(1,2,i,j),this%cell%A(2,2,i,j)], [2,2])
        this%cell%iA(:,:,i,j) = matmul(G, A) / det(this%cell%A(:,:,i,j))
        end associate
      end do
    end do
    ! Vertex
    do j = this%half_jms, this%half_jme
      do i = this%half_ims, this%half_ime
        associate (p0 => this%vtx  %p(:,i  ,j  ), &
                   p1 => this%yedge%p(:,i-1,j  ), &
                   p2 => this%yedge%p(:,i  ,j  ), &
                   p3 => this%xedge%p(:,i  ,j-1), &
                   p4 => this%xedge%p(:,i  ,j  ))
        e1 = cross_product(cross_product(p1, p2), p0); e1 = e1 / norm2(e1)
        e2 = cross_product(cross_product(p3, p4), p0); e2 = e2 / norm2(e2)
        this%vtx%sina(i,j) = dot_product(p0, cross_product(e1, e2))
        this%vtx%cosa(i,j) = dot_product(e1, e2)
        end associate
      end do
    end do
    ! X edge
    do j = this%full_jms, this%full_jme
      do i = this%half_ims, this%half_ime
        associate (p0 => this%xedge%p(:,i  ,j  ), &
                   p1 => this%cell %p(:,i-1,j  ), &
                   p2 => this%cell %p(:,i  ,j  ), &
                   p3 => this%vtx  %p(:,i  ,j  ), &
                   p4 => this%vtx  %p(:,i  ,j+1))
        e1 = cross_product(cross_product(p1, p2), p0); e1 = e1 / norm2(e1)
        e2 = cross_product(cross_product(p3, p4), p0); e2 = e2 / norm2(e2)
        this%xedge%sina(i,j) = dot_product(p0, cross_product(e1, e2))
        this%xedge%cosa(i,j) = dot_product(e1, e2)
        sinlon = sin(this%xedge%lon(i,j)); coslon = cos(this%xedge%lon(i,j))
        sinlat = sin(this%xedge%lat(i,j)); coslat = cos(this%xedge%lat(i,j))
        G  = reshape([1.0_r8,  this%xedge%cosa(i,j),  this%xedge%cosa(i,j), 1.0_r8], [2,2])
        iG = reshape([1.0_r8, -this%xedge%cosa(i,j), -this%xedge%cosa(i,j), 1.0_r8], [2,2]) / this%xedge%sina(i,j)**2
        elon = [-sinlon, coslon, 0.0_r8]
        elat = [-sinlat*coslon, -sinlat*sinlon, coslat]
        A = reshape([dot_product(e1,elon),dot_product(e1,elat),dot_product(e2,elon),dot_product(e2,elat)], [2,2])
        this%xedge%A(:,:,i,j) = matmul(A, iG)
        A = reshape([this%xedge%A(1,1,i,j),-this%xedge%A(2,1,i,j),-this%xedge%A(1,2,i,j),this%xedge%A(2,2,i,j)], [2,2])
        this%xedge%iA(:,:,i,j) = matmul(G, A) / det(this%xedge%A(:,:,i,j))
        end associate
      end do
    end do
    ! Y edge
    do j = this%half_jms, this%half_jme
      do i = this%full_ims, this%full_ime
        associate (p0 => this%yedge%p(:,i  ,j  ), &
                   p1 => this%vtx  %p(:,i  ,j  ), &
                   p2 => this%vtx  %p(:,i+1,j  ), &
                   p3 => this%cell %p(:,i  ,j-1), &
                   p4 => this%cell %p(:,i  ,j  ))
        e1 = cross_product(cross_product(p1, p2), p0); e1 = e1 / norm2(e1)
        e2 = cross_product(cross_product(p3, p4), p0); e2 = e2 / norm2(e2)
        this%yedge%sina(i,j) = dot_product(p0, cross_product(e1, e2))
        this%yedge%cosa(i,j) = dot_product(e1, e2)
        sinlon = sin(this%yedge%lon(i,j)); coslon = cos(this%yedge%lon(i,j))
        sinlat = sin(this%yedge%lat(i,j)); coslat = cos(this%yedge%lat(i,j))
        G  = reshape([1.0_r8,  this%yedge%cosa(i,j),  this%yedge%cosa(i,j), 1.0_r8], [2,2])
        iG = reshape([1.0_r8, -this%yedge%cosa(i,j), -this%yedge%cosa(i,j), 1.0_r8], [2,2]) / this%yedge%sina(i,j)**2
        elon = [-sinlon, coslon, 0.0_r8]
        elat = [-sinlat*coslon, -sinlat*sinlon, coslat]
        A = reshape([dot_product(e1,elon),dot_product(e1,elat),dot_product(e2,elon),dot_product(e2,elat)], [2,2])
        this%yedge%A(:,:,i,j) = matmul(A, iG)
        A = reshape([this%yedge%A(1,1,i,j),-this%yedge%A(2,1,i,j),-this%yedge%A(1,2,i,j),this%yedge%A(2,2,i,j)], [2,2])
        this%yedge%iA(:,:,i,j) = matmul(G, A) / det(this%yedge%A(:,:,i,j))
        end associate
      end do
    end do

    ! --------------------------------------------------------------------------
    ! Setup areas.
    !
    !    p4____p34______p3
    !     |  4  |  3  |
    !  p41|_____|p0___|p23
    !     |  1  |  2  |
    !     |_____|_____|
    !    p1    p12    p2
    !
    ! Cell
    do j = this%full_jms, this%full_jme
      do i = this%full_ims, this%full_ime
        associate (p0  => this%cell %p(:,i  ,j  ), &
                   p1  => this%vtx  %p(:,i  ,j  ), &
                   p2  => this%vtx  %p(:,i+1,j  ), &
                   p3  => this%vtx  %p(:,i+1,j+1), &
                   p4  => this%vtx  %p(:,i  ,j+1), &
                   p12 => this%yedge%p(:,i  ,j  ), &
                   p23 => this%xedge%p(:,i+1,j  ), &
                   p34 => this%yedge%p(:,i  ,j+1), &
                   p41 => this%xedge%p(:,i  ,j))
        this%cell%area   (  i,j) = spherical_rectangle_area(radius, p1 , p2 , p3 , p4 )
        this%cell%subarea(1,i,j) = spherical_rectangle_area(radius, p1 , p12, p0 , p41)
        this%cell%subarea(2,i,j) = spherical_rectangle_area(radius, p12, p2 , p23, p0 )
        this%cell%subarea(3,i,j) = spherical_rectangle_area(radius, p0 , p23, p3 , p34)
        this%cell%subarea(4,i,j) = spherical_rectangle_area(radius, p41, p0 , p34, p4 )
        end associate
      end do
    end do
    ! Vertex
    do j = this%half_jms, this%half_jme
      do i = this%half_ims, this%half_ime
        associate (p0  => this%vtx  %p(:,i  ,j  ), &
                   p1  => this%cell %p(:,i-1,j-1), &
                   p2  => this%cell %p(:,i  ,j-1), &
                   p3  => this%cell %p(:,i  ,j  ), &
                   p4  => this%cell %p(:,i-1,j  ), &
                   p12 => this%xedge%p(:,i  ,j-1), &
                   p23 => this%yedge%p(:,i  ,j  ), &
                   p34 => this%xedge%p(:,i  ,j  ), &
                   p41 => this%yedge%p(:,i-1,j  ))
        this%vtx%area   (  i,j) = spherical_rectangle_area(radius, p1 , p2 , p3 , p4 )
        this%vtx%subarea(1,i,j) = spherical_rectangle_area(radius, p1 , p12, p0 , p41)
        this%vtx%subarea(2,i,j) = spherical_rectangle_area(radius, p12, p2 , p23, p0 )
        this%vtx%subarea(3,i,j) = spherical_rectangle_area(radius, p0 , p23, p3 , p34)
        this%vtx%subarea(4,i,j) = spherical_rectangle_area(radius, p41, p0 , p34, p4 )
        end associate
      end do
    end do
    !           p3
    !           /|\
    !          / | \
    !         /  |  \
    !        /   |   \
    !     p4 \   |   /p2
    !         \  |  /
    !          \ | /
    !           \|/
    !           p1
    !
    ! X edge
    do j = this%full_jms, this%full_jme
      do i = this%half_ims, this%half_ime
        associate (p1 => this%vtx %p(:,i  ,j  ), &
                   p2 => this%cell%p(:,i  ,j  ), &
                   p3 => this%vtx %p(:,i  ,j+1), &
                   p4 => this%cell%p(:,i-1,j  ))
        this%xedge%subarea(1,i,j) = spherical_triangle_area(radius, p1, p3, p4)
        this%xedge%subarea(2,i,j) = spherical_triangle_area(radius, p1, p2, p3)
        this%xedge%area(i,j) = sum(this%xedge%subarea(1:2,i,j))
        end associate
      end do
    end do
    !           p3
    !           /\
    !          /  \
    !         /    \
    !      p4/______\ p2
    !        \      /
    !         \    /
    !          \  /
    !           \/
    !           p1
    ! Y edge
    do j = this%half_jms, this%half_jme
      do i = this%full_ims, this%full_ime
        associate (p1 => this%cell%p(:,i  ,j-1), &
                   p2 => this%vtx %p(:,i+1,j  ), &
                   p3 => this%cell%p(:,i  ,j  ), &
                   p4 => this%vtx %p(:,i  ,j  ))
        this%yedge%subarea(1,i,j) = spherical_triangle_area(radius, p1, p2, p4)
        this%yedge%subarea(2,i,j) = spherical_triangle_area(radius, p2, p3, p4)
        this%yedge%area(i,j) = sum(this%yedge%subarea(1:2,i,j))
        end associate
      end do
    end do

    ! Setup neighbor panels.
    if (this%active) then
      even_panel = mod(this%id, 2) == 0
      if (this%full_ids == 1) then ! Left neighbor
        ngb_ids = merge(nc - hw + 1      , nc - this%full_jde + 1, even_panel)
        ngb_jds = merge(this%full_jds    , nc - hw + 1           , even_panel)
        ngb_ndx = merge(hw               , this%full_ndy         , even_panel)
        ngb_ndy = merge(this%full_ndy    , hw                    , even_panel)
        call this%ngb(1)%ptr%init(proj_type=proj_type, nc=nc, nz=nz, hw=hw, ids=ngb_ids, jds=ngb_jds, ndx=ngb_ndx, ndy=ngb_ndy)
      end if
      if (this%full_ide == this%nx) then ! Right neighbor
        ngb_ids = merge(nc - this%full_jde + 1, 1                , even_panel)
        ngb_jds = merge(1                     , this%full_jds    , even_panel)
        ngb_ndx = merge(this%full_ndy         , hw               , even_panel)
        ngb_ndy = merge(hw                    , this%full_ndy    , even_panel)
        call this%ngb(2)%ptr%init(proj_type=proj_type, nc=nc, nz=nz, hw=hw, ids=ngb_ids, jds=ngb_jds, ndx=ngb_ndx, ndy=ngb_ndy)
      end if
      if (this%full_jds == 1) then ! Bottom neighbor
        ngb_ids = merge(nc - hw + 1           , this%full_ids    , even_panel)
        ngb_jds = merge(nc - this%full_ide + 1, nc - hw + 1      , even_panel)
        ngb_ndx = merge(hw                    , this%full_ndx    , even_panel)
        ngb_ndy = merge(this%full_ndx         , hw               , even_panel)
        call this%ngb(3)%ptr%init(proj_type=proj_type, nc=nc, nz=nz, hw=hw, ids=ngb_ids, jds=ngb_jds, ndx=ngb_ndx, ndy=ngb_ndy)
      end if
      if (this%full_jde == this%ny) then ! Top neighbor
        ngb_ids = merge(this%full_ids    , 1                     , even_panel)
        ngb_jds = merge(1                , nc - this%full_ide + 1, even_panel)
        ngb_ndx = merge(this%full_ndx    , hw                    , even_panel)
        ngb_ndy = merge(hw               , this%full_ndx         , even_panel)
        call this%ngb(4)%ptr%init(proj_type=proj_type, nc=nc, nz=nz, hw=hw, ids=ngb_ids, jds=ngb_jds, ndx=ngb_ndx, ndy=ngb_ndy)
      end if
    end if

    this%initialized = .true.

  end subroutine cubed_sphere_panel_init

  subroutine cubed_sphere_panel_clear(this)

    class(cubed_sphere_panel_type), intent(inout) :: this

    call this%cell %clear()
    call this%vtx  %clear()
    call this%xedge%clear()
    call this%yedge%clear()

    this%initialized = .false.

  end subroutine cubed_sphere_panel_clear

  subroutine cubed_sphere_panel_connect(this, id, left_panel, right_panel, bottom_panel, top_panel)

    class(cubed_sphere_panel_type), intent(inout) :: this
    integer, intent(in) :: id ! Panel ID (1-6)
    type(cubed_sphere_panel_type), intent(in), target :: left_panel
    type(cubed_sphere_panel_type), intent(in), target :: right_panel
    type(cubed_sphere_panel_type), intent(in), target :: bottom_panel
    type(cubed_sphere_panel_type), intent(in), target :: top_panel

    this%id = id
    this%ngb(1)%ptr => left_panel
    this%ngb(2)%ptr => right_panel
    this%ngb(3)%ptr => bottom_panel
    this%ngb(4)%ptr => top_panel

  end subroutine cubed_sphere_panel_connect

  subroutine cubed_sphere_panel_final(this)

    type(cubed_sphere_panel_type), intent(inout) :: this

    call this%clear()

  end subroutine cubed_sphere_panel_final

  subroutine grid_init(this, is, ie, js, je)

    class(grid_type), intent(inout) :: this
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: js
    integer, intent(in) :: je

    call this%clear()

    allocate(this%p      (  3,is:ie,js:je)); this%p         = inf
    allocate(this%x1     (    is:ie      )); this%x1        = inf
    allocate(this%y1     (          js:je)); this%y1        = inf
    allocate(this%lx     (    is:ie,js:je)); this%lx        = inf
    allocate(this%ly     (    is:ie,js:je)); this%ly        = inf
    allocate(this%lon    (    is:ie,js:je)); this%lon       = inf
    allocate(this%lon_deg(    is:ie,js:je)); this%lon_deg   = inf
    allocate(this%lat    (    is:ie,js:je)); this%lat       = inf
    allocate(this%lat_deg(    is:ie,js:je)); this%lat_deg   = inf
    allocate(this%area   (    is:ie,js:je)); this%area      = inf
    allocate(this%subarea(  4,is:ie,js:je)); this%subarea   = inf
    allocate(this%sina   (    is:ie,js:je)); this%sina      = inf
    allocate(this%cosa   (    is:ie,js:je)); this%cosa      = inf
    allocate(this%A      (2,2,is:ie,js:je)); this%A         = inf
    allocate(this%iA     (2,2,is:ie,js:je)); this%iA        = inf

  end subroutine grid_init

  subroutine grid_clear(this)

    class(grid_type), intent(inout) :: this

    if (allocated(this%p      )) deallocate(this%p      )
    if (allocated(this%x1     )) deallocate(this%x1     )
    if (allocated(this%y1     )) deallocate(this%y1     )
    if (allocated(this%lx     )) deallocate(this%lx     )
    if (allocated(this%ly     )) deallocate(this%ly     )
    if (allocated(this%lon    )) deallocate(this%lon    )
    if (allocated(this%lon_deg)) deallocate(this%lon_deg)
    if (allocated(this%lat    )) deallocate(this%lat    )
    if (allocated(this%lat_deg)) deallocate(this%lat_deg)
    if (allocated(this%area   )) deallocate(this%area   )
    if (allocated(this%subarea)) deallocate(this%subarea)
    if (allocated(this%sina   )) deallocate(this%sina   )
    if (allocated(this%cosa   )) deallocate(this%cosa   )
    if (allocated(this%A      )) deallocate(this%A      )
    if (allocated(this%iA     )) deallocate(this%iA     )

  end subroutine grid_clear

  subroutine grid_set_lon_lat(this, ipn)

    class(grid_type), intent(inout) :: this
    integer, intent(in) :: ipn

    real(r8), parameter :: one = 1.0_r8
    integer i, j

    do j = lbound(this%y1, 1), ubound(this%y1, 1)
      do i = lbound(this%x1, 1), ubound(this%x1, 1)
        associate (x1 => this%x1(i), y1 => this%y1(j), p => this%p(:,i,j), lon => this%lon(i,j), lat => this%lat(i,j))
        select case (ipn)
        case (1)
          p = [ one,   x1,   y1]
        case (2)
          p = [ -x1,  one,   y1]
        case (3)
          p = [ -x1,  -y1,  one]
        case (4)
          p = [-one,  -y1,  -x1]
        case (5)
          p = [  y1, -one,  -x1]
        case (6)
          p = [  y1,   x1, -one]
        end select
        p = p / norm2(p)
        call xyz2lonlat(1.0d0, p(1), p(2), p(3), lon, lat)
        this%lon_deg(i,j) = lon * deg
        this%lat_deg(i,j) = lat * deg
      end associate
      end do
    end do

  end subroutine grid_set_lon_lat

  subroutine grid_final(this)

    type(grid_type), intent(inout) :: this

    call this%clear()

  end subroutine grid_final

end module cubed_sphere_panel_mod
