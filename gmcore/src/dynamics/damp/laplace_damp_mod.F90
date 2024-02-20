module laplace_damp_mod

  ! ∂q       n+1    2n
  ! -- = (-1)    K ∇ q
  ! ∂t

  use const_mod
  use latlon_mesh_mod, only: global_mesh
  use latlon_field_types_mod
  use latlon_parallel_mod
  use process_mod, only: proc

  implicit none

  private

  public laplace_damp_init
  public laplace_damp_final
  public laplace_damp_run

  interface laplace_damp_run
    module procedure laplace_damp_run_2d
  end interface laplace_damp_run

contains

  subroutine laplace_damp_init()

  end subroutine laplace_damp_init

  subroutine laplace_damp_final()

  end subroutine laplace_damp_final

  subroutine laplace_damp_run_2d(f, order, coef)

    type(latlon_field2d_type), intent(inout) :: f
    integer, intent(in) :: order
    real(r8), intent(in), optional :: coef

    type(latlon_field2d_type) g1, g2, fx, fy
    real(r8) work(f%mesh%full_ids:f%mesh%full_ide), pole
    real(r8) c0, s
    integer i, j, k, l

    c0 = 1.0_r8; if (present(coef)) c0 = coef
    s = (-1)**(order / 2)

    select case (f%loc)
    case ('cell')
      call g1%init('', '', '', 'cell', f%mesh, f%halo)
      call g2%init('', '', '', 'cell', f%mesh, f%halo)
      call fx%init('', '', '', 'lon' , f%mesh, f%halo)
      call fy%init('', '', '', 'lat' , f%mesh, f%halo)
      g1%d = f%d
      do l = 1, (order - 2) / 2
        do j = f%mesh%full_jds_no_pole, f%mesh%full_jde_no_pole
          do i = f%mesh%half_ids - 1, f%mesh%half_ide
            fx%d(i,j) = (g1%d(i+1,j) - g1%d(i,j)) / f%mesh%de_lon(j)
          end do
        end do
        do j = f%mesh%half_jds - merge(0, 1, f%mesh%has_south_pole()), f%mesh%half_jde
          do i = f%mesh%full_ids, f%mesh%full_ide
            fy%d(i,j) = (g1%d(i,j+1) - g1%d(i,j)) / f%mesh%de_lat(j)
          end do
        end do
        do j = f%mesh%full_jds_no_pole, f%mesh%full_jde_no_pole
          do i = f%mesh%full_ids, f%mesh%full_ide
            g2%d(i,j) = (                                    &
              (fx%d(i,j) - fx%d(i-1,j)) * f%mesh%le_lon(j) + &
               fy%d(i,j  ) * f%mesh%le_lat(j  ) -            &
               fy%d(i,j-1) * f%mesh%le_lat(j-1)              &
            ) / f%mesh%area_cell(j)
          end do
        end do
        if (f%mesh%has_south_pole()) then
          j = f%mesh%full_jds
          do i = f%mesh%full_ids, f%mesh%full_ide
            work(i) = fy%d(i,j)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * f%mesh%le_lat(j) / global_mesh%full_nlon / f%mesh%area_cell(j)
          do i = f%mesh%full_ids, f%mesh%full_ide
            g2%d(i,j) = pole
          end do
        end if
        if (f%mesh%has_north_pole()) then
          j = f%mesh%full_jde
          do i = f%mesh%full_ids, f%mesh%full_ide
            work(i) = fy%d(i,j-1)
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = -pole * f%mesh%le_lat(j-1) / global_mesh%full_nlon / f%mesh%area_cell(j)
          do i = f%mesh%full_ids, f%mesh%full_ide
            g2%d(i,j) = pole
          end do
        end if
        call fill_halo(g2)
        g1%d = g2%d
      end do
      do j = f%mesh%full_jds_no_pole, f%mesh%full_jde_no_pole
        do i = f%mesh%half_ids - 1, f%mesh%half_ide
          fx%d(i,j) = (g1%d(i+1,j) - g1%d(i,j)) / f%mesh%de_lon(j)
        end do
      end do
      do j = f%mesh%half_jds - merge(0, 1, f%mesh%has_south_pole()), f%mesh%half_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          fy%d(i,j) = (g1%d(i,j+1) - g1%d(i,j)) / f%mesh%de_lat(j)
        end do
      end do
      if (order > 2) then
        do j = f%mesh%full_jds_no_pole, f%mesh%full_jde_no_pole
          do i = f%mesh%half_ids - 1, f%mesh%half_ide
            fx%d(i,j) = fx%d(i,j) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j) * (f%d(i+1,j) - f%d(i,j))))
          end do
        end do
        do j = f%mesh%half_jds - merge(0, 1, f%mesh%has_south_pole()), f%mesh%half_jde
          do i = f%mesh%full_ids, f%mesh%full_ide
            fy%d(i,j) = fy%d(i,j) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j) * (f%d(i,j+1) - f%d(i,j))))
          end do
        end do
      end if
      do j = f%mesh%full_jds_no_pole, f%mesh%full_jde_no_pole
        do i = f%mesh%full_ids, f%mesh%full_ide
          f%d(i,j) = f%d(i,j) - s * c0 * (       &
            (                                    &
              fx%d(i,j) - fx%d(i-1,j)            &
            ) * f%mesh%le_lon(j) +               &
            (                                    &
              fy%d(i,j  ) * f%mesh%le_lat(j  ) - &
              fy%d(i,j-1) * f%mesh%le_lat(j-1)   &
            )) / f%mesh%area_cell(j)
        end do
      end do
      if (f%mesh%has_south_pole()) then
        j = f%mesh%full_jds
        do i = f%mesh%full_ids, f%mesh%full_ide
          work(i) = fy%d(i,j)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = s * c0 * pole * f%mesh%le_lat(j) / global_mesh%full_nlon / f%mesh%area_cell(j)
        do i = f%mesh%full_ids, f%mesh%full_ide
          f%d(i,j) = f%d(i,j) - pole
        end do
      end if
      if (f%mesh%has_north_pole()) then
        j = f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          work(i) = fy%d(i,j-1)
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -s * c0 * pole * f%mesh%le_lat(j-1) / global_mesh%full_nlon / f%mesh%area_cell(j)
        do i = f%mesh%full_ids, f%mesh%full_ide
          f%d(i,j) = f%d(i,j) - pole
        end do
      end if
    case ('lev')

    end select
    call g1%clear()
    call g2%clear()
    call fx%clear()
    call fy%clear()

  end subroutine laplace_damp_run_2d

end module laplace_damp_mod
