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
!   This module setups Gaussian convolution filter for each zonal circle to
!   improve the polar stability of lat-lon grid.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module filter_types_mod

  use const_mod
  use namelist_mod
  use math_mod
  use latlon_mesh_mod

  implicit none

  private

  public filter_type

  type filter_type
    type(latlon_mesh_type), pointer :: mesh => null()
    real(8), allocatable :: width_lon(:)
    integer, allocatable :: ngrid_lon(:)
    real(8), allocatable :: wgt_lon(:,:)
    real(8), allocatable :: width_lat(:)
    integer, allocatable :: ngrid_lat(:)
    real(8), allocatable :: wgt_lat(:,:)
  contains
    procedure :: init  => filter_init
    procedure :: clear => filter_clear
    final :: filter_final
  end type filter_type

contains

  subroutine gaussian_weight(width, ngrid, w)

    real(8), intent(in) :: width
    integer, intent(in) :: ngrid
    real(8), intent(out) :: w(:)

    real(8) s
    integer i, x

    s = width / filter_gauss_sigma
    do i = 1, ngrid
      x = i - (ngrid + 1) / 2
      w(i) = exp(-x**2 / (2 * s**2)) / (s * sqrt(pi2))
    end do
    w(:ngrid) = w(:ngrid) / sum(w(:ngrid))

  end subroutine gaussian_weight

  subroutine filter_init(this, mesh, type)

    class(filter_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: type

    real(8) dx, dt, cfl, w
    integer j, n

    call this%clear()

    this%mesh => mesh
    dt = dt_dyn
    allocate(this%width_lon(mesh%full_jms:mesh%full_jme)); this%width_lon = 0
    allocate(this%ngrid_lon(mesh%full_jms:mesh%full_jme)); this%ngrid_lon = 0
    allocate(this%width_lat(mesh%half_jms:mesh%half_jme)); this%width_lat = 0
    allocate(this%ngrid_lat(mesh%half_jms:mesh%half_jme)); this%ngrid_lat = 0

    if (max_wave_speed > 0) then
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        dx = mesh%de_lon(j)
        if (dx > 0) then
          cfl = max_wave_speed * dt / dx
          w = filter_coef_a * cfl / max_cfl * (filter_coef_b * (tanh(90 - abs(mesh%full_lat_deg(j))) - 1) + 1)
          w = max(filter_min_width, w)
          n = ceiling(w); if (mod(n, 2) == 0) n = n + 1; n = max(3, n)
          this%width_lon(j) = w
          this%ngrid_lon(j) = n
        end if
      end do
      do j = mesh%half_jds, mesh%half_jde
        dx = merge(mesh%de_lon(j+1), mesh%de_lon(j), mesh%half_lat(j) < 0)
        if (dx > 0) then
          cfl = max_wave_speed * dt / dx
          w = filter_coef_a * cfl / max_cfl * (filter_coef_b * (tanh(90 - abs(mesh%half_lat_deg(j))) - 1) + 1)
          w = max(filter_min_width, w)
          n = ceiling(w); if (mod(n, 2) == 0) n = n + 1; n = max(3, n)
          this%width_lat(j) = w
          this%ngrid_lat(j) = n
        end if
      end do
    end if

    allocate(this%wgt_lon(maxval(this%ngrid_lon),mesh%full_jms:mesh%full_jme)); this%wgt_lon = 0
    allocate(this%wgt_lat(maxval(this%ngrid_lat),mesh%half_jms:mesh%half_jme)); this%wgt_lat = 0
    select case (type)
    case ('big_filter')
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        if (this%ngrid_lon(j) > 1) then
          call gaussian_weight(this%width_lon(j), this%ngrid_lon(j), this%wgt_lon(:,j))
        end if
      end do
      do j = mesh%half_jds, mesh%half_jde
        if (this%ngrid_lat(j) > 1) then
          call gaussian_weight(this%width_lat(j), this%ngrid_lat(j), this%wgt_lat(:,j))
        end if
      end do
    case ('small_filter')
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        if (this%ngrid_lon(j) > 1) then
          w = filter_coef_c * this%width_lon(j)
          n = ceiling(w); if (mod(n, 2) == 0) n = n + 1; n = max(3, n)
          this%width_lon(j) = w
          this%ngrid_lon(j) = min(n, this%ngrid_lon(j))
          call gaussian_weight(this%width_lon(j), this%ngrid_lon(j), this%wgt_lon(:,j))
        end if
      end do
      do j = mesh%half_jds, mesh%half_jde
        if (this%ngrid_lat(j) > 1) then
          w = filter_coef_c * this%width_lat(j)
          n = ceiling(w); if (mod(n, 2) == 0) n = n + 1; n = max(3, n)
          this%width_lat(j) = w
          this%ngrid_lat(j) = min(n, this%ngrid_lat(j))
          call gaussian_weight(this%width_lat(j), this%ngrid_lat(j), this%wgt_lat(:,j))
        end if
      end do
    case default
      call log_error('Invalid filter type ' // trim(type) // '!')
    end select

  end subroutine filter_init

  subroutine filter_clear(this)

    class(filter_type), intent(inout) :: this

    if (allocated(this%width_lon)) deallocate(this%width_lon)
    if (allocated(this%ngrid_lon)) deallocate(this%ngrid_lon)
    if (allocated(this%wgt_lon  )) deallocate(this%wgt_lon  )
    if (allocated(this%width_lat)) deallocate(this%width_lat)
    if (allocated(this%ngrid_lat)) deallocate(this%ngrid_lat)
    if (allocated(this%wgt_lat  )) deallocate(this%wgt_lat  )

  end subroutine filter_clear

  subroutine filter_final(this)

    type(filter_type), intent(inout) :: this

    call this%clear()

  end subroutine filter_final

end module filter_types_mod
