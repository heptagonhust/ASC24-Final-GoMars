module laplace_damp_mod

  ! ∂q       n+1    2n
  ! -- = (-1)    K ∇ q
  ! ∂t

  use const_mod
  use namelist_mod
  use math_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public laplace_damp_init
  public laplace_damp_final
  public laplace_damp_on_cell
  public laplace_damp_on_lon_edge
  public laplace_damp_on_lat_edge
  public decay_from_pole
  public decay_from_sfc
  public decay_from_top

  interface laplace_damp_on_cell
    module procedure laplace_damp_on_cell_2d
    module procedure laplace_damp_on_cell_3d
  end interface laplace_damp_on_cell

  integer, parameter :: diff_hw(1:8) = [1, 1, 2, 2, 3, 3, 4, 4]

  real(r8), target :: diff_weights(9,1:8) = reshape([  &
    -1,  1,   0,   0,   0,   0,   0,   0,  0,  & ! 1
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 8])

  real(r8), allocatable, dimension(:), target :: lat_ones
  real(r8), allocatable, dimension(:), target :: lev_ones
  real(r8), allocatable, dimension(:), target :: decay_from_pole
  real(r8), allocatable, dimension(:), target :: decay_from_sfc
  real(r8), allocatable, dimension(:), target :: decay_from_top

contains

  subroutine laplace_damp_init()

    integer k, k0

    allocate(lat_ones       (global_mesh%full_nlat)); lat_ones  = 1
    allocate(lev_ones       (global_mesh%full_nlev)); lev_ones  = 1
    allocate(decay_from_pole(global_mesh%full_nlat)); decay_from_pole = 1
    allocate(decay_from_sfc (global_mesh%full_nlev))
    allocate(decay_from_top (global_mesh%full_nlev))

    k0 = 5
    do k = global_mesh%full_kds, global_mesh%full_kde
      decay_from_sfc(k) = exp((global_mesh%full_nlev - k)**2 * log(0.1d0) / k0**2)
    end do

    do k = global_mesh%full_kds, global_mesh%full_kde
      decay_from_top(k) = exp_two_values(1.0_r8, 0.0_r8, 1.0_r8, 10.0_r8, real(k, r8))
    end do

  end subroutine laplace_damp_init

  subroutine laplace_damp_final()

    if (allocated(lat_ones       )) deallocate(lat_ones       )
    if (allocated(lev_ones       )) deallocate(lev_ones       )
    if (allocated(decay_from_pole)) deallocate(decay_from_pole)
    if (allocated(decay_from_sfc )) deallocate(decay_from_sfc )

  end subroutine laplace_damp_final

  subroutine laplace_damp_on_cell_2d(mesh, halo, order, f, coef, lat_coef, fill)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    logical, intent(in), optional :: fill

    real(r8) gx1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gx2(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gy1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gy2(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fx (mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fy (mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) c0, s, cj_half
    real(r8), pointer :: cj(:)
    integer i, j, k

    c0 = 0.5_r8**order * merge(coef, 1.0_r8, present(coef))
    if (present(lat_coef)) then
      cj => lat_coef
    else
      cj => lat_ones
    end if
    s = (-1)**(order / 2)

    gx1 = f
    gy1 = f
    do k = 1, (order - 2) / 2
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gx2(i,j) = gx1(i-1,j) - 2 * gx1(i,j) + gx1(i+1,j)
          gy2(i,j) = gy1(i,j-1) - 2 * gy1(i,j) + gy1(i,j+1)
        end do
      end do
      call fill_halo(halo, gx2, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
      call fill_halo(halo, gy2, full_lon=.true., full_lat=.true., west_halo=.false., east_halo=.false.)
      gx1 = gx2
      gy1 = gy2
    end do
    ! Calculate damping flux at interfaces.
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids - 1, mesh%half_ide
        fx(i,j) = s * cj(j) * (gx1(i+1,j) - gx1(i,j))
      end do
    end do
    do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
      cj_half = merge(cj(j), cj(j+1), mesh%half_lat(j) < 0)
      do i = mesh%full_ids, mesh%full_ide
        fy(i,j) = s * cj_half * (gy1(i,j+1) - gy1(i,j))
      end do
    end do
    ! Limit damping flux to avoid upgradient (Xue 2000).
    if (order > 2) then
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          fx(i,j) = fx(i,j) * max(0.0_r8, sign(1.0_r8, -fx(i,j) * (f(i+1,j) - f(i,j))))
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          fy(i,j) = fy(i,j) * max(0.0_r8, sign(1.0_r8, -fy(i,j) * (f(i,j+1) - f(i,j))))
        end do
      end do
    end if
    ! Update variable.
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%full_ids, mesh%full_ide
        f(i,j) = f(i,j) - c0 * (fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))
      end do
    end do
    if (merge(fill, .true., present(fill))) call fill_halo(halo, f, full_lon=.true., full_lat=.true.)

  end subroutine laplace_damp_on_cell_2d

  subroutine laplace_damp_on_cell_3d(mesh, halo, order, f, coef, lat_coef, lev_coef)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional :: lat_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional :: lev_coef(global_mesh%full_nlev)

    real(r8) ck
    integer k

    ck = 1
    do k = mesh%full_kds, mesh%full_kde
      if (present(coef) .and. present(lev_coef)) then
        ck = coef * lev_coef(k)
      else if (present(coef)) then
        ck = coef
      else if (present(lev_coef)) then
        ck = lev_coef(k)
      end if
      call laplace_damp_on_cell_2d(mesh, halo, order, f(:,:,k), coef=ck, lat_coef=lat_coef, fill=.false.)
    end do
    call fill_halo(halo, f, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine laplace_damp_on_cell_3d

  subroutine laplace_damp_on_lon_edge(mesh, halo, order, f, coef, lat_coef, lev_coef)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%half_ims:mesh%half_ime, &
                                 mesh%full_jms:mesh%full_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lev_coef(global_mesh%full_nlev)

    real(r8) gx1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gy1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gx2(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) gy2(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fx (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fy (mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme)
    real(r8) c0, cj_half, s
    real(r8), pointer :: cj(:), ck(:)
    integer i, j, k, l

    c0 = 0.5_r8**order * merge(coef, 1.0_r8, present(coef))
    if (present(lat_coef)) then
      cj => lat_coef
    else
      cj => lat_ones
    end if
    if (present(lev_coef)) then
      ck => lev_coef
    else
      ck => lev_ones
    end if
    s = (-1)**(order / 2)

    do k = mesh%full_kds, mesh%full_kde
      gx1 = f(:,:,k)
      gy1 = f(:,:,k)
      do l = 1, (order - 2) / 2
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            gx2(i,j) = gx1(i-1,j) - 2 * gx1(i,j) + gx1(i+1,j)
            gy2(i,j) = gy1(i,j-1) - 2 * gy1(i,j) + gy1(i,j+1)
          end do
        end do
        call fill_halo(halo, gx2, full_lon=.false., full_lat=.true., south_halo=.false., north_halo=.false.)
        call fill_halo(halo, gy2, full_lon=.false., full_lat=.true., west_halo=.false., east_halo=.false.)
        gx1 = gx2
        gy1 = gy2
      end do
      ! Calculate damping flux at interfaces.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide + 1
          fx(i,j) = s * cj(j) * (gx1(i,j) - gx1(i-1,j))
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        cj_half = merge(cj(j), cj(j+1), mesh%half_lat(j) < 0)
        do i = mesh%half_ids, mesh%half_ide
          fy(i,j) = s * cj_half * (gy1(i,j+1) - gy1(i,j))
        end do
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      if (order > 2) then
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide + 1
            fx(i,j) = fx(i,j) * max(0.0_r8, sign(1.0_r8, -fx(i,j) * (f(i,j,k) - f(i-1,j,k))))
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%half_ids, mesh%half_ide
            fy(i,j) = fy(i,j) * max(0.0_r8, sign(1.0_r8, -fy(i,j) * (f(i,j+1,k) - f(i,j,k))))
          end do
        end do
      end if
      ! Update variable.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          f(i,j,k) = f(i,j,k) - c0 * ck(k) * (fx(i+1,j) - fx(i,j) + fy(i,j) - fy(i,j-1))
        end do
      end do
    end do
    call fill_halo(halo, f, full_lon=.false., full_lat=.true., full_lev=.true.)

  end subroutine laplace_damp_on_lon_edge

  subroutine laplace_damp_on_lat_edge(mesh, halo, order, f, coef, lat_coef, lev_coef)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%half_jms:mesh%half_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lev_coef(global_mesh%full_nlev)

    real(r8) gx1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) gy1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) gx2(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) gy2(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) fx (mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme)
    real(r8) fy (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) c0, cj_half, s
    real(r8), pointer :: cj(:), ck(:)
    integer i, j, k, l

    c0 = 0.5_r8**order * merge(coef, 1.0_r8, present(coef))
    if (present(lat_coef)) then
      cj => lat_coef
    else
      cj => lat_ones
    end if
    if (present(lev_coef)) then
      ck => lev_coef
    else
      ck => lev_ones
    end if
    s = (-1)**(order / 2)

    do k = mesh%full_kds, mesh%full_kde
      gx1 = f(:,:,k)
      gy1 = f(:,:,k)
      do l = 1, (order - 2) / 2
        if (mesh%has_south_pole()) then
          j = mesh%half_jds
          gy1(:,j-1) = gy1(:,j)
        end if
        if (mesh%has_north_pole()) then
          j = mesh%half_jde
          gy1(:,j+1) = gy1(:,j)
        end if
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            gx2(i,j) = gx1(i-1,j) - 2 * gx1(i,j) + gx1(i+1,j)
            gy2(i,j) = gy1(i,j-1) - 2 * gy1(i,j) + gy1(i,j+1)
          end do
        end do
        call fill_halo(halo, gx2, full_lon=.true., full_lat=.false., south_halo=.false., north_halo=.false.)
        call fill_halo(halo, gy2, full_lon=.true., full_lat=.false., west_halo=.false., east_halo=.false.)
        gx1 = gx2
        gy1 = gy2
      end do
      if (mesh%has_south_pole()) then
        j = mesh%half_jds
        gy1(:,j-1) = gy1(:,j)
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_jde
        gy1(:,j+1) = gy1(:,j)
      end if
      ! Calculate damping flux at interfaces.
      do j = mesh%half_jds, mesh%half_jde
        cj_half = merge(cj(j), cj(j+1), mesh%half_lat(j) < 0)
        do i = mesh%half_ids - 1, mesh%half_ide
          fx(i,j) = s * cj_half * (gx1(i+1,j) - gx1(i,j))
        end do
      end do
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide
          fy(i,j) = s * cj(j) * (gy1(i,j) - gy1(i,j-1))
        end do
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      if (order > 2) then
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids - 1, mesh%half_ide
            fx(i,j) = fx(i,j) * max(0.0_r8, sign(1.0_r8, -fx(i,j) * (f(i+1,j,k) - f(i,j,k))))
          end do
        end do
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            fy(i,j) = fy(i,j) * max(0.0_r8, sign(1.0_r8, -fy(i,j) * (f(i,j,k) - f(i,j-1,k))))
          end do
        end do
      end if
      ! Update variable.
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          f(i,j,k) = f(i,j,k) - c0 * ck(k) * (fx(i,j) - fx(i-1,j) + fy(i,j+1) - fy(i,j))
        end do
      end do
    end do
    call fill_halo(halo, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine laplace_damp_on_lat_edge

end module laplace_damp_mod
