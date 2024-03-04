module latlon_operators_mod

  use const_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

contains

  subroutine div_operator(fx, fy, div, with_halo)

    type(latlon_field3d_type), intent(in   ) :: fx
    type(latlon_field3d_type), intent(in   ) :: fy
    type(latlon_field3d_type), intent(inout) :: div
    logical, intent(in), optional :: with_halo

    real(r8) work(div%mesh%full_ids:div%mesh%full_ide,div%nlev)
    real(r8) pole(div%nlev)
    integer i, j, k, is, ie, js, je, ks, ke

    associate (mesh => div%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, div%loc == 'cell')
    ke = merge(mesh%full_kde, mesh%half_kde, div%loc == 'cell')
    is = mesh%full_ids
    ie = mesh%full_ide; if (present(with_halo)) ie = merge(ie + 1, ie, with_halo)
    js = mesh%full_jds_no_pole
    je = mesh%full_jde_no_pole; if (present(with_halo)) je = merge(je + 1, je, with_halo)
    do k = ks, ke
      do j = js, je
        do i = is, ie
          div%d(i,j,k) = ((                    &
            fx%d(i,j,k) - fx%d(i-1,j,k)        &
          ) * mesh%le_lon(j) + (               &
            fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
            fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = fy%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = ks, ke
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = -pole(k)
        end do
      end do
    end if
    end associate

  end subroutine div_operator

end module latlon_operators_mod
