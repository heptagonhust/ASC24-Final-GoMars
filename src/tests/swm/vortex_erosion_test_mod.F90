module vortex_erosion_test_mod
  
  ! J. R. Bates & Yong Li (1997) Simulation of Stratospheric Vortex Erosion 
  ! Using Three Different Global Shallow Water Numerical Models.

  use flogger
  use string
  use const_mod
  use time_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public vortex_erosion_test_set_ic
  public vortex_erosion_test_apply_forcing

contains
  
  subroutine vortex_erosion_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, neval, ierr
    real(8) gh0, gz_, abserr
    
    associate(mesh => block%mesh           , &
              u    => block%dstate(1)%u_lon, &
              v    => block%dstate(1)%v_lat, &
              gz   => block%dstate(1)%gz   )
    gh0 = g * 6.0d3

    do j = mesh%full_jds, mesh%full_jde
      u(:,j,1) = u_function(mesh%full_lat(j))
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true.)
    
    v = 0

    do j = mesh%full_jds, mesh%full_jde
      i = mesh%half_ids
      if (j == mesh%full_jds) then
        gz(i,j,1) = gh0
      else
        call qags(gh_integrand, -0.5d0*pi, mesh%full_lat(j), 1.0d-12, 1.0d-3, gz_, abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_str(i) // ',' // to_str(j) // ')!', __FILE__, __LINE__)
        end if 
        gz(i,j,1) = gh0 - gz_
      end if
      do i = mesh%full_ids, mesh%full_ide
        gz(i,j,1) = gz(mesh%half_ids,j,1)
      end do 
    end do 

    call fill_halo(block%halo, gz, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine vortex_erosion_test_set_ic

  subroutine vortex_erosion_test_apply_forcing(block, static)

    type(block_type), intent(in) :: block
    type(static_type), intent(inout) :: static
    integer i, j, k
    real(8) hs, elapsed_days, at, b_lat, y
    
    hs = 720.0d0 ! m
    elapsed_days = elapsed_seconds / 86400.0d0
    elapsed_days = mod(elapsed_days, 20.0d0)

    associate (mesh => block%mesh, &
               gzs  => static%gzs)
    if (elapsed_days < 4) then
      at = 0.5d0 * (1 - cos(pi * elapsed_days / 4.0d0))
    else if (elapsed_days < 16) then
      at = 1.0d0
    else if (elapsed_days < 20) then
      at = 0.5d0 * (1 + cos(pi * (elapsed_days - 16) / 4.0d0))
    end if

    b_lat = 0
    gzs = 0
    do j = mesh%full_jds, mesh%full_jde
      if (mesh%full_lat(j) > 0) then
        y = (tan(pi * 0.25d0) / tan(mesh%full_lat(j)))**2
        b_lat = y * exp(1 - y)
        do i = mesh%full_ids, mesh%full_ide
          gzs(i,j) = hs * at * b_lat * mesh%full_sin_lon(i) * g
        end do
      end if
    end do
    call fill_halo(block%halo, gzs, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine vortex_erosion_test_apply_forcing

  real(8) function gh_integrand(lat) result(res)

    real(8), intent(in) :: lat

    real(8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(8) function u_function(lat) result(res)

    real(8), intent(in) :: lat
    real(8) lat_deg

    lat_deg = lat * deg

    if (lat_deg <= 0) then
      res = -lat_deg / 9.0d0 - 10
    else if (lat_deg < 60) then
      res = lat_deg - 10
    else
      res = -5.0d0 / 3.0d0 * lat_deg + 150
    end if 

  end function u_function

end module vortex_erosion_test_mod
