! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module interp_mod

  use const_mod
  use namelist_mod
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use latlon_field_types_mod
  use perf_mod

  implicit none

  private

  !                              / lev_lat_edge
  !               o-------------o------------o lev_vtx
  !              /|            /            /|
  !             / |                        / |
  !            /  |        |              /  |
  !           o   |        o lev_edge   -o- lev_lon_edge
  !          /    |        |            /    |
  !         /     o vtx                /     o vtx
  !        /      |                   /      |
  !       o-------+-----o------------o       |
  !       |       |                  |       |
  ! lon_edge -o-  |        o cell    |  -o- lon_edge
  !       |       |                  |       |
  !       |       o------------------+-------o
  !       |      /       /           |      /
  !       o vtx /       o lat_edge   o vtx /
  !       |    /       /             |    /
  !       |   o                      |   o
  !       |  /                       |  /
  !       | /                        | /
  !       |/                         |/
  !       o-------------o------------o
  !

  public interp_init
  public interp_final
  public interp_run
  public average_run
  public interp_cell_to_height_level
  public interp_lon_edge_to_height_level
  public interp_lat_edge_to_height_level
  public interp_lev_edge_to_height_level
  public interp_cell_to_pressure_level
  public interp_lon_edge_to_pressure_level
  public interp_lat_edge_to_pressure_level
  public interp_lev_edge_to_pressure_level

  interface interp_run
    module procedure interp_run_3d
  end interface interp_run

  interface average_run
    module procedure average_run_3d
  end interface average_run

contains

  subroutine interp_init()

    call interp_final()

  end subroutine interp_init

  subroutine interp_final()

  end subroutine interp_final

  subroutine interp_run_3d(x, y)

    type(latlon_field3d_type), intent(in   ) :: x
    type(latlon_field3d_type), intent(inout) :: y

    real(r8) a, b, c, x1, x2, x3
    integer i, j, k

    select case (trim(x%loc) // '>' // trim(y%loc))
    ! --------------------------------------------------------------------------
    case ('cell>lon')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%mesh%area_lon_west(j) * x%d(i  ,j,k) + &
                          x%mesh%area_lon_east(j) * x%d(i+1,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lat')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_north(j) * x%d(i,j+1,k) + &
                          x%mesh%area_lat_south(j) * x%d(i,j  ,k)   &
                         ) / x%mesh%area_lat(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lev')
      if (x%mesh%full_nlev == 1) return
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! -------
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (2 * x%mesh%half_dlev(k))
        b = x%mesh%full_dlev(k  ) / (2 * x%mesh%half_dlev(k))
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k)
          end do
        end do
      end do
      k = x%mesh%half_kds
      x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
      x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1)
        end do
      end do
      k = x%mesh%half_kde
      x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
      x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2)
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>vtx')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (                                                   &
              (x%d(i,j  ,k) + x%d(i+1,j  ,k)) * x%mesh%area_subcell(2,j  ) + &
              (x%d(i,j+1,k) + x%d(i+1,j+1,k)) * x%mesh%area_subcell(1,j+1)   &
            ) / x%mesh%area_vtx(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lon>cell')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lon_east(j) * x%d(i-1,j,k) + &
                          x%mesh%area_lon_west(j) * x%d(i  ,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lat>cell')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_south(j  ) * x%d(i,j  ,k) + &
                          x%mesh%area_lat_north(j-1) * x%d(i,j-1,k)   &
                         ) / (x%mesh%area_lat_south(j) + x%mesh%area_lat_north(j-1))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>cell')
      ! =======
      !
      ! ---o--- k
      !
      ! ===?=== k
      !
      ! ---o--- k+1
      !
      ! =======
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = 0.5_r8 * (x%d(i,j,k) + x%d(i,j,k+1))
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lon>lev_lon')
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! ----o--
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        do j = x%mesh%full_jds, x%mesh%full_jde
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k-1)
          end do
        end do
      end do
      k = x%mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
      x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
      x3 = x%mesh%full_lev(k+2) - x%mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%half_ids, x%mesh%half_ide
          y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1) + c * x%d(i,j,k+2)
        end do
      end do
      k = x%mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
      x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
      x3 = x%mesh%half_lev(k) - x%mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%full_jds, x%mesh%full_jde
        do i = x%mesh%half_ids, x%mesh%half_ide
          y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2) + c * x%d(i,j,k-3)
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lat>lev_lat')
      ! -------
      !
      ! ===o=== k-1
      !
      ! ---?--- k
      !
      ! ===o=== k
      !
      ! -------
      do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
        a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k-1)
          end do
        end do
      end do
      k = x%mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
      x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
      x3 = x%mesh%full_lev(k+2) - x%mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%half_jds, x%mesh%half_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k) + b * x%d(i,j,k+1) + c * x%d(i,j,k+2)
        end do
      end do
      k = x%mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = x%mesh%half_lev(k) - x%mesh%full_lev(k-1)
      x2 = x%mesh%half_lev(k) - x%mesh%full_lev(k-2)
      x3 = x%mesh%half_lev(k) - x%mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = x%mesh%half_jds, x%mesh%half_jde
        do i = x%mesh%full_ids, x%mesh%full_ide
          y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2) + c * x%d(i,j,k-3)
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>lev_lon')
      do k = x%mesh%half_kds, x%mesh%half_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%mesh%area_lon_west(j) * x%d(i  ,j,k) + &
                          x%mesh%area_lon_east(j) * x%d(i+1,j,k)   &
                         ) / x%mesh%area_lon(j)
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('lev>lev_lat')
      do k = x%mesh%half_kds, x%mesh%half_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%mesh%area_lat_north(j) * x%d(i,j+1,k) + &
                          x%mesh%area_lat_south(j) * x%d(i,j  ,k)   &
                         ) / x%mesh%area_lat(j)
          end do
        end do
      end do
    end select

  end subroutine interp_run_3d

  subroutine average_run_3d(x, y)

    type(latlon_field3d_type), intent(in   ) :: x
    type(latlon_field3d_type), intent(inout) :: y

    integer i, j, k

    ! real(r8) tempx(x%mesh%full_jds:x%mesh%full_jde+1, x%mesh%full_ids:x%mesh%full_ide, x%mesh%full_kds:x%mesh%full_kde)

    call t_startf ('average_run')

    select case (trim(x%loc) // '>' // trim(y%loc))
    ! --------------------------------------------------------------------------
    case ('cell>lon')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
          do i = x%mesh%half_ids, x%mesh%half_ide
            y%d(i,j,k) = (x%d(i,j,k) + x%d(i+1,j,k)) * 0.5_r8
          end do
        end do
      end do
    ! --------------------------------------------------------------------------
    case ('cell>lat')
      do k = x%mesh%full_kds, x%mesh%full_kde
        do j = x%mesh%half_jds, x%mesh%half_jde
          do i = x%mesh%full_ids, x%mesh%full_ide
            y%d(i,j,k) = (x%d(i,j,k) + x%d(i,j+1,k)) * 0.5_r8
          end do
        end do
      end do
    !   do k = x%mesh%full_kds, x%mesh%full_kde
    !     do j = x%mesh%half_jds, x%mesh%half_jde+1
    !       do i = x%mesh%full_ids, x%mesh%full_ide
    !         tempx(j,i,k) = x%d(i,j,k)
    !       end do
    !     end do
    !   end do

    !   do k = x%mesh%full_kds, x%mesh%full_kde
    !       do i = x%mesh%full_ids, x%mesh%full_ide
    !         do j = x%mesh%half_jds, x%mesh%half_jde
    !         y%d(i,j,k) = (tempx(j,i,k) + tempx(j+1,i,k)) * 0.5_r8
    !       end do
    !     end do
    !   end do
    end select

    call t_stopf ('average_run')

  end subroutine average_run_3d

  subroutine interp_cell_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: zo
    real(r8), intent(inout) :: y(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme)

    real(r8) dz1, dz2, z1, z2, a, b
    integer i, j, k

    ! --o-- z(k-1)
    !
    ! --?-- zo
    !
    ! --o-- z(k)
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            y(i,j) = (dz2 * x(i,j,k-1) + dz1 * x(i,j,k)) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_kde) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k) + b * x(i,j,k-1)
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_kds + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k-1) + b * x(i,j,k)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_cell_to_height_level

  subroutine interp_lon_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%half_ims:mesh%half_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    !               x1
    ! x(i-1,k-1) o--x--o x(i,k-1)
    !
    !               ? zo
    !
    ! x(i-1,k  ) o--x--o x(i,k  )
    !               x2
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = (dz2 * x1 + dz1 * x2) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_kde) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_kds + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_height_level

  subroutine interp_lat_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%half_jms:mesh%half_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = (dz2 * x1 + dz1 * x2) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_kde) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_kds + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_height_level

  subroutine interp_lev_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%half_kms:mesh%half_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%half_kms:mesh%half_kme)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%half_kde, mesh%half_kds + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            y(i,j) = (dz2 * x(i,j,k-1) + dz1 * x(i,j,k)) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_kde) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k) + b * x(i,j,k-1)
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_kds + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k-1) + b * x(i,j,k)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_height_level

  subroutine interp_cell_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: po
    real(r8), intent(inout) :: y(mesh%full_ims:mesh%full_ime, &
                                 mesh%full_jms:mesh%full_jme)
    logical, intent(in), optional :: logp

    logical logp_opt
    real(r8) p0, dp1, dp2
    integer i, j, k

    logp_opt = .false.; if (present(logp)) logp_opt = logp

    p0 = merge(log(po), po, logp_opt)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_opt) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            y(i,j) = (dp2 * x(i,j,k-1) + dp1 * x(i,j,k)) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_cell_to_pressure_level

  subroutine interp_lon_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%half_ims:mesh%half_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)
    logical, intent(in), optional :: logp

    logical logp_opt
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_opt = .false.; if (present(logp)) logp_opt = logp

    p0 = merge(log(po), po, logp_opt)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_opt) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = (dp2 * x1 + dp1 * x2) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_pressure_level

  subroutine interp_lat_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%half_jms:mesh%half_jme, &
                              mesh%full_kms:mesh%full_kme)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)
    logical, intent(in), optional :: logp

    logical logp_opt
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_opt = .false.; if (present(logp)) logp_opt = logp

    p0 = merge(log(po), po, logp_opt)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%full_kde, mesh%full_kds + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_opt) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = (dp2 * x1 + dp1 * x2) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_pressure_level

  subroutine interp_lev_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%half_kms:mesh%half_kme)
    real(r8), intent(in) :: x(mesh%full_ims:mesh%full_ime, &
                              mesh%full_jms:mesh%full_jme, &
                              mesh%half_kms:mesh%half_kme)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme)
    logical, intent(in), optional :: logp

    logical logp_opt
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_opt = .false.; if (present(logp)) logp_opt = logp

    p0 = merge(log(po), po, logp_opt)

    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%half_kde, mesh%half_kds + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_opt) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            y(i,j) = (dp2 * x(i,j,k-1) + dp1 * x(i,j,k)) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_pressure_level

  subroutine interp_lon_edge_to_cell(mesh, x_lon, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lon(mesh%half_ims:mesh%half_ime, &
                                  mesh%full_jms:mesh%full_jme, &
                                  mesh%full_kms:mesh%full_kme)
    real(r8), intent(out) :: x(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme, &
                               mesh%full_kms:mesh%full_kme)

    integer i, j, k

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          x(i,j,k) = (mesh%area_lon_east(j) * x_lon(i-1,j,k) + &
                      mesh%area_lon_west(j) * x_lon(i  ,j,k)   &
                     ) / mesh%area_lon(j)
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_cell

  subroutine interp_lat_edge_to_cell(mesh, x_lat, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lat(mesh%full_ims:mesh%full_ime, &
                                  mesh%half_jms:mesh%half_jme, &
                                  mesh%full_kms:mesh%full_kme)
    real(r8), intent(out) :: x(mesh%full_ims:mesh%full_ime, &
                               mesh%full_jms:mesh%full_jme, &
                               mesh%full_kms:mesh%full_kme)

    integer i, j, k

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          x(i,j,k) = (mesh%area_lat_south(j  ) * x_lat(i,j  ,k) + &
                      mesh%area_lat_north(j-1) * x_lat(i,j-1,k)   &
                     ) / (mesh%area_lat_south(j) + mesh%area_lat_north(j-1))
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_cell

  subroutine interp_lon_edge_to_lev_lon_edge(mesh, x_lon, x_lev_lon, handle_top_bottom)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lon(mesh%half_ims:mesh%half_ime, &
                                  mesh%full_jms:mesh%full_jme, &
                                  mesh%full_kms:mesh%full_kme)
    real(r8), intent(out) :: x_lev_lon(mesh%half_ims:mesh%half_ime, &
                                       mesh%full_jms:mesh%full_jme, &
                                       mesh%half_kms:mesh%half_kme)
    logical, intent(in), optional :: handle_top_bottom

    logical handle_top_bottom_opt
    integer i, j, k
    real(r8) x1, x2, x3, a, b, c

    handle_top_bottom_opt = .false.; if (present(handle_top_bottom)) handle_top_bottom_opt = handle_top_bottom

    ! -------
    !
    ! ===o=== k-1
    !
    ! ---?--- k
    !
    ! ===o=== k
    !
    ! ----o--
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      a = mesh%full_dlev(k-1) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      b = mesh%full_dlev(k  ) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          x_lev_lon(i,j,k) = a * x_lon(i,j,k) + b * x_lon(i,j,k-1)
        end do
      end do
    end do

    if (handle_top_bottom_opt) then
      k = mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      x3 = mesh%full_lev(k+2) - mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          x_lev_lon(i,j,k) = a * x_lon(i,j,k) + b * x_lon(i,j,k+1) + c * x_lon(i,j,k+2)
        end do
      end do
      k = mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      x3 = mesh%half_lev(k) - mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%half_ids, mesh%half_ide
          x_lev_lon(i,j,k) = a * x_lon(i,j,k-1) + b * x_lon(i,j,k-2) + c * x_lon(i,j,k-3)
        end do
      end do
    end if

  end subroutine interp_lon_edge_to_lev_lon_edge

  subroutine interp_lat_edge_to_lev_lat_edge(mesh, x_lat, x_lev_lat, handle_top_bottom)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lat(mesh%full_ims:mesh%full_ime, &
                                  mesh%half_jms:mesh%half_jme, &
                                  mesh%full_kms:mesh%full_kme)
    real(r8), intent(out) :: x_lev_lat(mesh%full_ims:mesh%full_ime, &
                                       mesh%half_jms:mesh%half_jme, &
                                       mesh%half_kms:mesh%half_kme)
    logical, intent(in), optional :: handle_top_bottom

    logical handle_top_bottom_opt
    integer i, j, k
    real(r8) x1, x2, x3, a, b, c

    handle_top_bottom_opt = .false.; if (present(handle_top_bottom)) handle_top_bottom_opt = handle_top_bottom

    ! -------
    !
    ! ===o=== k-1
    !
    ! ---?--- k
    !
    ! ===o=== k
    !
    ! -------
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      a = mesh%full_dlev(k-1) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      b = mesh%full_dlev(k  ) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          x_lev_lat(i,j,k) = a * x_lat(i,j,k) + b * x_lat(i,j,k-1)
        end do
      end do
    end do

    if (handle_top_bottom_opt) then
      k = mesh%half_kds
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      x3 = mesh%full_lev(k+2) - mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          x_lev_lat(i,j,k) = a * x_lat(i,j,k) + b * x_lat(i,j,k+1) + c * x_lat(i,j,k+2)
        end do
      end do
      k = mesh%half_kde
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      x3 = mesh%half_lev(k) - mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          x_lev_lat(i,j,k) = a * x_lat(i,j,k-1) + b * x_lat(i,j,k-2) + c * x_lat(i,j,k-3)
        end do
      end do
    end if

  end subroutine interp_lat_edge_to_lev_lat_edge

end module interp_mod
