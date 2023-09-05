module laplace_damp_mod

  ! ∂q       n+1    2n
  ! -- = (-1)    K ∇ q
  ! ∂t

  use const_mod
  use namelist_mod
  use math_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use block_mod

  implicit none

  private

  public laplace_damp_init
  public laplace_damp_final
  public laplace_damp_on_cell
  public laplace_damp_on_lon_edge
  public laplace_damp_on_lat_edge
  public lev_ones
  public decay_from_pole
  public decay_from_sfc
  public decay_from_top

  interface laplace_damp_on_cell
    module procedure laplace_damp_on_cell_2d
    module procedure laplace_damp_on_cell_3d
  end interface laplace_damp_on_cell

  real(r8), allocatable, dimension(:), target :: lat_ones
  real(r8), allocatable, dimension(:), target :: lev_ones
  real(r8), allocatable, dimension(:), target :: decay_from_pole
  real(r8), allocatable, dimension(:), target :: decay_from_sfc
  real(r8), allocatable, dimension(:), target :: decay_from_top

contains

  subroutine laplace_damp_init()

    integer j, j0, k, k0

    allocate(lat_ones       (global_mesh%full_nlat)); lat_ones  = 1
    allocate(lev_ones       (global_mesh%full_nlev)); lev_ones  = 1
    allocate(decay_from_pole(global_mesh%full_nlat)); decay_from_pole = 0
    allocate(decay_from_sfc (global_mesh%full_nlev))
    allocate(decay_from_top (global_mesh%full_nlev))

    do j = global_mesh%full_jds_no_pole, global_mesh%full_jde_no_pole
      decay_from_pole(j) = exp_two_values(0.01_r8, 1.0_r8, real(abs(global_mesh%full_lat_deg(2)), r8), 50.0_r8, real(abs(global_mesh%full_lat_deg(j)), r8))
    end do

    k0 = 3
    do k = global_mesh%full_kds, global_mesh%full_kde
      decay_from_sfc(k) = exp_two_values(1.0_r8, 0.0_r8, real(global_mesh%full_kde, r8), real(k0, r8), real(k, r8))
    end do

    k0 = 20
    do k = global_mesh%full_kds, global_mesh%full_kde
      decay_from_top(k) = exp_two_values(1.0_r8, 0.0_r8, 1.0_r8, real(k0, r8), real(k, r8))
    end do

  end subroutine laplace_damp_init

  subroutine laplace_damp_final()

    if (allocated(lat_ones       )) deallocate(lat_ones       )
    if (allocated(lev_ones       )) deallocate(lev_ones       )
    if (allocated(decay_from_pole)) deallocate(decay_from_pole)
    if (allocated(decay_from_sfc )) deallocate(decay_from_sfc )

  end subroutine laplace_damp_final

  subroutine laplace_damp_on_cell_2d(mesh, halo, order, f, coef, lon_coef, lat_coef, fill)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lon_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    logical, intent(in), optional :: fill

    real(r8) g1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) g2(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fx(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fy(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) c0, s
    real(r8) work(mesh%full_ids:mesh%full_ide), pole
    real(r8), pointer :: ci(:), cj(:)
    logical fill_opt
    integer i, j, l

    c0 = 1.0_r8; if (present(coef)) c0 = coef
    if (present(lon_coef)) then
      ci => lon_coef
    else
      ci => lat_ones
    end if
    if (present(lat_coef)) then
      cj => lat_coef
    else
      cj => lat_ones
    end if
    s = (-1)**(order / 2)

    g1 = f
    do l = 1, (order - 2) / 2
      ! Here we consider 2nd-order diffusion or damping.
      ! Firstly, calculate damping flux which is a gradient operator.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          fx(i,j) = (g1(i+1,j) - g1(i,j)) / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          fy(i,j) = (g1(i,j+1) - g1(i,j)) / mesh%de_lat(j)
        end do
      end do
      ! Secondly, calculate 2nd-order damping which is a divergence operator on damping flux.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          g2(i,j) = (                                &
            (fx(i,j) - fx(i-1,j)) * mesh%le_lon(j) + &
             fy(i,j  ) * mesh%le_lat(j  ) -          &
             fy(i,j-1) * mesh%le_lat(j-1)            &
          ) / mesh%area_cell(j)
        end do
      end do
      ! Handle poles.
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do i = mesh%full_ids, mesh%full_ide
          work(i) = fy(i,j)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do i = mesh%full_ids, mesh%full_ide
          g2(i,j) = pole
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          work(i) = fy(i,j-1)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do i = mesh%full_ids, mesh%full_ide
          g2(i,j) = pole
        end do
      end if
      call fill_halo(halo, g2, full_lon=.true., full_lat=.true.)
      ! Copy values for next iteration.
      g1 = g2
    end do
    ! Calculate final damping flux at interfaces.
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids - 1, mesh%half_ide
        fx(i,j) = (g1(i+1,j) - g1(i,j)) / mesh%de_lon(j)
      end do
    end do
    do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        fy(i,j) = (g1(i,j+1) - g1(i,j)) / mesh%de_lat(j)
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
        f(i,j) = f(i,j) - s * c0 * (       &
          ci(j) * (                        &
            fx(i,j) - fx(i-1,j)            &
          ) * mesh%le_lon(j) +             &
          cj(j) * (                        &
            fy(i,j  ) * mesh%le_lat(j  ) - &
            fy(i,j-1) * mesh%le_lat(j-1)   &
          )                                &
        ) / mesh%area_cell(j)
      end do
    end do
    ! Handle poles.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do i = mesh%full_ids, mesh%full_ide
        work(i) = fy(i,j)
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = s * c0 * cj(j) * pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do i = mesh%full_ids, mesh%full_ide
        f(i,j) = f(i,j) - pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        work(i) = fy(i,j-1)
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = -s * c0 * cj(j) * pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do i = mesh%full_ids, mesh%full_ide
        f(i,j) = f(i,j) - pole
      end do
    end if
    fill_opt = .true.; if (present(fill)) fill_opt = fill
    if (fill_opt) call fill_halo(halo, f, full_lon=.true., full_lat=.true.)

  end subroutine laplace_damp_on_cell_2d

  subroutine laplace_damp_on_cell_3d(mesh, halo, order, f, coef, lon_coef, lat_coef, lev_coef)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional :: lon_coef(global_mesh%full_nlat)
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
      call laplace_damp_on_cell_2d(mesh, halo, order, f(:,:,k), coef=ck, lon_coef=lon_coef, lat_coef=lat_coef, fill=.false.)
    end do
    call fill_halo(halo, f, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)

  end subroutine laplace_damp_on_cell_3d

  subroutine laplace_damp_on_lon_edge(mesh, halo, order, f, coef, lon_coef, lat_coef, lev_coef, fill)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%half_ims:mesh%half_ime, &
                                 mesh%full_jms:mesh%full_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lon_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lev_coef(global_mesh%full_nlev)
    logical, intent(in), optional :: fill

    real(r8) g1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) g2(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fx(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) fy(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme)
    real(r8) c0, s
    real(r8) work(mesh%half_ids:mesh%half_ide), pole
    real(r8), pointer :: ci(:), cj(:), ck(:)
    logical fill_opt
    integer i, j, k, l

    c0 = 1.0_r8; if (present(coef)) c0 = coef
    if (present(lon_coef)) then
      ci => lon_coef
    else
      ci => lat_ones
    end if
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
      g1 = f(:,:,k)
      do l = 1, (order - 2) / 2
        ! Here we consider 2nd-order diffusion or damping.
        ! Firstly, calculate damping flux which is a gradient operator.
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide + 1
            fx(i,j) = (g1(i,j) - g1(i-1,j)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy(i,j) = (g1(i,j+1) - g1(i,j)) / mesh%de_lat(j)
          end do
        end do
        ! Secondly, calculate 2nd-order damping which is a divergence operator on damping flux.
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            g2(i,j) = (                                &
              (fx(i+1,j) - fx(i,j)) * mesh%le_lon(j) + &
               fy(i,j  ) * mesh%le_lat(j  ) -          &
               fy(i,j-1) * mesh%le_lat(j-1)            &
            ) / mesh%area_lon(j) / 2
          end do
        end do
        ! Handle poles.
        if (mesh%has_south_pole()) then
          j = mesh%full_jds
          do i = mesh%half_ids, mesh%half_ide
            work(i) = fy(i,j)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
          do i = mesh%half_ids, mesh%half_ide
            g2(i,j) = pole
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            work(i) = fy(i,j-1)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = -pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
          do i = mesh%half_ids, mesh%half_ide
            g2(i,j) = pole
          end do
        end if
        call fill_halo(halo, g2, full_lon=.false., full_lat=.true.)
        ! Copy values for next iteration.
        g1 = g2
      end do
      ! Calculate damping flux at interfaces.
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide + 1
          fx(i,j) = (g1(i,j) - g1(i-1,j)) / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          fy(i,j) = (g1(i,j+1) - g1(i,j)) / mesh%de_lat(j)
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
          f(i,j,k) = f(i,j,k) - s * c0 * (   &
            ci(j) * (                        &
              fx(i+1,j) - fx(i,j)            &
            ) * mesh%le_lon(j) +             &
            cj(j) * (                        &
              fy(i,j  ) * mesh%le_lat(j  ) - &
              fy(i,j-1) * mesh%le_lat(j-1)   &
            )                                &
          ) / mesh%area_lon(j) / 2
        end do
      end do
    end do
    fill_opt = .true.; if (present(fill)) fill_opt = fill
    if (fill_opt) call fill_halo(halo, f, full_lon=.false., full_lat=.true., full_lev=.true.)

  end subroutine laplace_damp_on_lon_edge

  subroutine laplace_damp_on_lat_edge(mesh, halo, order, f, coef, lon_coef, lat_coef, lev_coef, fill)

    type(mesh_type), intent(in) :: mesh
    type(halo_type), intent(in) :: halo(:)
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(mesh%full_ims:mesh%full_ime, &
                                 mesh%half_jms:mesh%half_jme, &
                                 mesh%full_kms:mesh%full_kme)
    real(r8), intent(in), optional :: coef
    real(r8), intent(in), optional, target :: lon_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lat_coef(global_mesh%full_nlat)
    real(r8), intent(in), optional, target :: lev_coef(global_mesh%full_nlev)
    logical, intent(in), optional :: fill

    real(r8) g1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) g2(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme)
    real(r8) fx(mesh%half_ims:mesh%half_ime,mesh%half_jms:mesh%half_jme)
    real(r8) fy(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme)
    real(r8) c0, ci_, cj_, s
    real(r8), pointer :: ci(:), cj(:), ck(:)
    logical fill_opt
    integer i, j, k, l

    c0 = 1.0_r8; if (present(coef)) c0 = coef
    if (present(lon_coef)) then
      ci => lon_coef
    else
      ci => lat_ones
    end if
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
      g1 = f(:,:,k)
      do l = 1, (order - 2) / 2
        ! Here we consider 2nd-order diffusion or damping.
        ! Firstly, calculate damping flux which is a gradient operator.
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids - 1, mesh%half_ide
            fx(i,j) = (g1(i+1,j) - g1(i,j)) / mesh%le_lat(j)
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            fy(i,j) = (g1(i,j) - g1(i,j-1)) / mesh%le_lon(j)
          end do
        end do
        ! Secondly, calculate 2nd-order damping which is a divergence operator on damping flux.
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            g2(i,j) = (                                &
              (fx(i,j) - fx(i-1,j)) * mesh%de_lat(j) + &
               fy(i,j+1) * mesh%de_lon(j+1) -          &
               fy(i,j  ) * mesh%de_lon(j  )            &
            ) / mesh%area_lat(j) / 2
          end do
        end do
        call fill_halo(halo, g2, full_lon=.true., full_lat=.false.)
        g1 = g2
      end do
      ! Calculate damping flux at interfaces.
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids - 1, mesh%half_ide
          fx(i,j) = (g1(i+1,j) - g1(i,j)) / mesh%le_lat(j)
        end do
      end do
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide
          fy(i,j) = (g1(i,j) - g1(i,j-1)) / mesh%le_lon(j)
        end do
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      if (order > 2) then
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids - 1, mesh%half_ide
            fx(i,j) = fx(i,j) * max(0.0_r8, sign(1.0_r8, -fx(i,j) * (f(i+1,j,k) - f(i,j,k))))
          end do
        end do
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            fy(i,j) = fy(i,j) * max(0.0_r8, sign(1.0_r8, -fy(i,j) * (f(i,j,k) - f(i,j-1,k))))
          end do
        end do
      end if
      ! Update variable.
      do j = mesh%half_jds, mesh%half_jde
        ci_ = merge(ci(j+1), ci(j), mesh%half_lat(j) < 0)
        cj_ = merge(cj(j+1), cj(j), mesh%half_lat(j) < 0)
        do i = mesh%full_ids, mesh%full_ide
          f(i,j,k) = f(i,j,k) - s * c0 * (   &
            ci_ * (                          &
              fx(i,j) - fx(i-1,j)            &
            ) * mesh%de_lat(j) +             &
            cj_ * (                          &
              fy(i,j+1) * mesh%le_lon(j+1) - &
              fy(i,j  ) * mesh%le_lon(j  )   &
            )                                &
          ) / mesh%area_lat(j) / 2
        end do
      end do
    end do
    fill_opt = .true.; if (present(fill)) fill_opt = fill
    if (fill_opt) call fill_halo(halo, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine laplace_damp_on_lat_edge

end module laplace_damp_mod
