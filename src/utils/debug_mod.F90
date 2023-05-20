module debug_mod

  use, intrinsic :: ieee_arithmetic
  use flogger
  use const_mod
  use namelist_mod
  use latlon_mesh_mod
  use block_mod
  use parallel_mod

  private

  public debug_check_areas
  public debug_print_min_max
  public debug_print_min_max_cell
  public debug_print_min_max_lev_edge
  public debug_print_min_max_lon_edge
  public debug_print_min_max_lat_edge
  public debug_is_inf

  interface debug_print_min_max
    module procedure debug_print_min_max_1d
    module procedure debug_print_min_max_2d
    module procedure debug_print_min_max_3d
  end interface debug_print_min_max

  interface debug_is_inf
    module procedure debug_is_inf_r4
    module procedure debug_is_inf_r8
  end interface debug_is_inf

contains

  subroutine debug_check_areas()

    type(latlon_mesh_type), pointer :: mesh
    real(8) total_area
    integer j

    mesh => global_mesh

    total_area = 0.0_r8
    do j = mesh%full_jds, mesh%full_jde
      total_area = total_area + mesh%area_cell(j) * mesh%full_nlon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate cell area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%half_jds, mesh%half_jde
      total_area = total_area + mesh%area_vtx(j) * mesh%half_nlon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate vertex area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%full_jds, mesh%full_jde
      total_area = total_area + sum(mesh%area_subcell(:,j)) * mesh%full_nlon * 2
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
    end if

    do j = mesh%full_jds, mesh%full_jde
      if (abs(mesh%area_cell(j) - 2.0d0 * sum(mesh%area_subcell(:,j))) / mesh%area_cell(j) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%half_jds, mesh%half_jde
      if (abs(mesh%area_vtx(j) - 2.0_r8 * (mesh%area_subcell(2,j) + mesh%area_subcell(1,j+1))) / mesh%area_vtx(j) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      if (abs(mesh%area_lon_north(j) + mesh%area_lon_south(j) - mesh%area_lon(j)) / mesh%area_lon(j) > 1.0d-12) then
        call log_error('Failed to calculate north and south subcell on lon grids!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%half_jds, mesh%half_jde
      if (abs(mesh%area_lat_west(j) + mesh%area_lat_east(j) - mesh%area_lat(j)) / mesh%area_lat(j) > 1.0d-11) then
        call log_error('Failed to calculate west and east subcell on lat grids!', __FILE__, __LINE__)
      end if
    end do

    total_area = 0.0_r8
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      total_area = total_area + mesh%area_lon(j) * mesh%full_nlon
    end do
    do j = mesh%half_jds, mesh%half_jde
      total_area = total_area + mesh%area_lat(j) * mesh%full_nlon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-9) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      total_area = total_area + (mesh%area_lon_north(j) + mesh%area_lon_south(j)) * mesh%full_nlon
    end do
    do j = mesh%half_jds, mesh%half_jde
      total_area = total_area + (mesh%area_lat_west(j) + mesh%area_lat_east(j)) * mesh%full_nlon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-9) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

  end subroutine debug_check_areas

  subroutine debug_print_min_max_1d(array, label)

    real(r8), intent(in) :: array(:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_1d

  subroutine debug_print_min_max_2d(array, label)

    real(r8), intent(in) :: array(:,:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_2d

  subroutine debug_print_min_max_3d(array, label)

    real(r8), intent(in) :: array(:,:,:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_3d

  subroutine debug_print_min_max_cell(mesh, array, label)

    type(latlon_mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_ims:mesh%full_ime, &
                                  mesh%full_jms:mesh%full_jme, &
                                  mesh%full_kms:mesh%full_kme)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%full_kds:mesh%full_kde)), &
                             maxval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%full_kds:mesh%full_kde))

  end subroutine debug_print_min_max_cell

  subroutine debug_print_min_max_lev_edge(mesh, array, label)

    type(latlon_mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_ims:mesh%full_ime, &
                                  mesh%full_jms:mesh%full_jme, &
                                  mesh%half_kms:mesh%half_kme)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%half_kds:mesh%half_kde)), &
                             maxval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%half_kds:mesh%half_kde))

  end subroutine debug_print_min_max_lev_edge

  subroutine debug_print_min_max_lon_edge(mesh, array, label)

    type(latlon_mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%half_ims:mesh%half_ime, &
                                  mesh%full_jms:mesh%full_jme, &
                                  mesh%full_kms:mesh%full_kme)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%half_ids:mesh%half_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%full_kds:mesh%full_kde)), &
                             maxval(array(mesh%half_ids:mesh%half_ide,   &
                                          mesh%full_jds:mesh%full_jde,   &
                                          mesh%full_kds:mesh%full_kde))

  end subroutine debug_print_min_max_lon_edge

  subroutine debug_print_min_max_lat_edge(mesh, array, label)

    type(latlon_mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_ims:mesh%full_ime, &
                                  mesh%half_jms:mesh%half_jme, &
                                  mesh%full_kms:mesh%full_kme)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%half_jds:mesh%half_jde,   &
                                          mesh%full_kds:mesh%full_kde)), &
                             maxval(array(mesh%full_ids:mesh%full_ide,   &
                                          mesh%half_jds:mesh%half_jde,   &
                                          mesh%full_kds:mesh%full_kde))

  end subroutine debug_print_min_max_lat_edge

  logical function debug_is_inf_r4(x) result(res)

    real(4), intent(in) :: x

    res = .not. ieee_is_finite(x)

  end function debug_is_inf_r4

  logical function debug_is_inf_r8(x) result(res)

    real(8), intent(in) :: x

    res = .not. ieee_is_finite(x)

  end function debug_is_inf_r8

end module debug_mod
