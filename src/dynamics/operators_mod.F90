module operators_mod

  use const_mod
  use vert_coord_mod
  use block_mod
  use parallel_mod
  use formula_mod
  use namelist_mod
  use log_mod
  use pgf_mod
  use adv_mod
  use nh_mod
  use interp_mod

  implicit none

  private

  public operators_init
  public operators_prepare
  public calc_ph
  public calc_m
  public calc_t
  public calc_gz_lev
  public calc_we_lev
  public calc_div
  public calc_vor
  public calc_coriolis
  public calc_grad_ke
  public calc_grad_mf
  public calc_grad_ptf
  public calc_dphsdt
  public calc_wedudlev_wedvdlev
  public nh_prepare
  public nh_solve

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

  interface
    subroutine interp_pv_interface(block, dstate, dt)
      import block_type, dstate_type, r8
      type(block_type), intent(inout) :: block
      type(dstate_type), intent(inout) :: dstate
      real(r8), intent(in) :: dt
    end subroutine interp_pv_interface
  end interface

  procedure(interp_pv_interface), pointer :: interp_pv => null()

contains

  subroutine operators_init()

    select case (pv_scheme)
    case ('midpoint')
      interp_pv => interp_pv_midpoint
    case ('upwind')
      interp_pv => interp_pv_upwind
    case ('tvd')
      interp_pv => interp_pv_tvd
    case default
      call log_error('Invalid pv_scheme ' // trim(pv_scheme) // '!', pid=proc%id)
    end select
    
  end subroutine operators_init

  subroutine operators_prepare_1(blocks, itime, dt)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    integer iblk

    do iblk = 1, size(blocks)
      if (baroclinic ) call calc_ph       (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_m                         (blocks(iblk), blocks(iblk)%dstate(itime))
      if (baroclinic ) call calc_t        (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_mf                        (blocks(iblk), blocks(iblk)%dstate(itime), dt)
      call calc_ke                        (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_pv                        (blocks(iblk), blocks(iblk)%dstate(itime))
      call interp_pv                      (blocks(iblk), blocks(iblk)%dstate(itime), dt)
      if (hydrostatic) call calc_gz_lev   (blocks(iblk), blocks(iblk)%dstate(itime))
      call pgf_prepare                    (blocks(iblk), blocks(iblk)%dstate(itime))
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, dstate, dt, pass)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    select case (pass)
    ! --------------------------------------------------------------------------
    case (all_pass)
      if (baroclinic ) call calc_ph       (block, dstate)
      call calc_m                         (block, dstate)
      if (baroclinic ) call calc_t        (block, dstate)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate)
      call calc_pv                        (block, dstate)
      call interp_pv                      (block, dstate, dt)
      if (hydrostatic) call calc_gz_lev   (block, dstate)
      call pgf_prepare                    (block, dstate)
    ! --------------------------------------------------------------------------
    case (forward_pass)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate)
      call calc_pv                        (block, dstate)
      call interp_pv                      (block, dstate, dt)
    ! --------------------------------------------------------------------------
    case (backward_pass)
      if (hydrostatic) then
        call calc_t                       (block, dstate)
        call calc_gz_lev                  (block, dstate)
      end if
      call pgf_prepare                    (block, dstate)
    end select

  end subroutine operators_prepare_2

  subroutine calc_ph(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    associate (mesh       => block%mesh       , &
               phs        => dstate%phs       , & ! in
               ph_lev     => dstate%ph_lev    , & ! out
               ph_exn_lev => dstate%ph_exn_lev, & ! out
               ph         => dstate%ph        )   ! out
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
          ph_exn_lev(i,j,k) = ph_lev(i,j,k)**Rd_o_cpd
        end do
      end do
    end do
    call fill_halo(block%halo, ph_exn_lev, full_lon=.true., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ph(i,j,k) = 0.5_r8 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
        end do
      end do
    end do
    call fill_halo(block%halo, ph, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine calc_ph

  subroutine calc_t(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    associate (mesh => block%mesh, &
               pt   => dstate%pt , & ! in
               ph   => dstate%ph , & ! in
               qv   => dstate%qv , & ! in
               t    => dstate%t  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          t(i,j,k) = temperature(pt(i,j,k), ph(i,j,k), qv(i,j,k))
        end do
      end do
    end do
    end associate

  end subroutine calc_t

  subroutine calc_we_lev(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(in) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k, l
    real(r8) mf

    associate (mesh       => block%mesh       , &
               dmfdlon    => dtend%dmfdlon    , & ! in
               dmfdlat    => dtend%dmfdlat    , & ! in
               dphs       => dtend%dphs       , & ! in
               m_lev      => dstate%m_lev     , & ! in
               we_lev     => dstate%we_lev    , & ! out
               we_lev_lon => dstate%we_lev_lon, & ! out
               we_lev_lat => dstate%we_lev_lat)   ! out
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mf = 0.0_r8
          do l = 1, k - 1
            mf = mf + dmfdlon(i,j,l) + dmfdlat(i,j,l)
          end do
          we_lev(i,j,k) = - vert_coord_calc_dphdt_lev(k, dphs(i,j)) - mf
        end do
      end do
    end do
    ! Set vertical boundary conditions.
    we_lev(:,:,mesh%half_kds) = 0.0_r8
    we_lev(:,:,mesh%half_kde) = 0.0_r8
    call fill_halo(block%halo, we_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    call block%adv_batch_pt%accum_we_lev(we_lev, m_lev, dt)

    call interp_lev_edge_to_lev_lon_edge(mesh, we_lev, we_lev_lon)
    call interp_lev_edge_to_lev_lat_edge(mesh, we_lev, we_lev_lat)
    end associate

  end subroutine calc_we_lev

  subroutine calc_ke(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) ke_vtx(4)
    real(r8) work(dstate%mesh%full_ids:dstate%mesh%full_ide,dstate%mesh%full_nlev)
    real(r8) pole(dstate%mesh%full_nlev)

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, & ! in
               v    => dstate%v_lat, & ! in
               ke   => dstate%ke   )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          ke(i,j,k) = (mesh%area_lon_west (j  ) * u(i-1,j  ,k)**2 + &
                       mesh%area_lon_east (j  ) * u(i  ,j  ,k)**2 + &
                       mesh%area_lat_north(j-1) * v(i  ,j-1,k)**2 + &
                       mesh%area_lat_south(j  ) * v(i  ,j  ,k)**2   &
                      ) / mesh%area_cell(j)
        end do
      end do
    end do

    if (ke_scheme == 2) then
      !
      !      ________u_________________u________
      !     |     i-1,j+1     |       i,j+1     |
      !     |                 |                 |
      !     |                 |                 |
      !     |        1        |        4        |
      !     v        o--------v--------o        v
      !  i-1,j    i-1,j      i,j      i,j    i+1,j
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     |________u________|________u________|
      !     |     i-1,j      i,j      i,j       |
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     v        o--------v--------o        v
      !  i-1,j-1  i-1,j-1    i,j-1    i,j-1  i+1,j-1
      !     |        2        |        3        |
      !     |                 |                 |
      !     |________u________|________u________|
      !           i-1,j-1             i,j-1
      !
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            ke_vtx(1) = (                                  &
              mesh%area_lat_east (j  ) * v(i-1,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v(i  ,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u(i-1,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u(i-1,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke_vtx(2) = (                                  &
              mesh%area_lat_east (j-1) * v(i-1,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v(i  ,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u(i-1,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u(i-1,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(3) = (                                  &
              mesh%area_lat_east (j-1) * v(i  ,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v(i+1,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u(i  ,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u(i  ,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(4) = (                                  &
              mesh%area_lat_east (j  ) * v(i  ,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v(i+1,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u(i  ,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u(i  ,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke(i,j,k) = (1.0_r8 - ke_cell_wgt) * (               &
              (ke_vtx(1) + ke_vtx(4)) * mesh%area_subcell(2,j) + &
              (ke_vtx(2) + ke_vtx(3)) * mesh%area_subcell(1,j)   &
            ) / mesh%area_cell(j) + ke_cell_wgt * ke(i,j,k)
          end do
        end do
      end do
    end if

    ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v(i,j,k)**2
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v(i,j-1,k)**2
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke(i,j,k) = pole(k)
        end do
      end do
    end if
    end associate

  end subroutine calc_ke

  subroutine calc_div(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    real(r8) work(dstate%mesh%full_ids:dstate%mesh%full_ide,dstate%mesh%full_nlev)
    real(r8) pole(dstate%mesh%full_nlev)
    integer i, j, k

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, & ! in
               v    => dstate%v_lat, & ! in
               div  => dstate%div  , & ! out
               div2 => dstate%div2 )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          div(i,j,k) = (                                                    &
            (u(i,j,k) * mesh%le_lon(j) - u(i-1,  j,k) * mesh%le_lon(j  )) + &
            (v(i,j,k) * mesh%le_lat(j) - v(i  ,j-1,k) * mesh%le_lat(j-1))   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          div(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -v(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          div(i,j,k) = pole(k)
        end do
      end do
    end if
    if (div_damp_order == 4) then
      call fill_halo(block%halo, div, full_lon=.true., full_lat=.true., full_lev=.true.)
    else
      call fill_halo(block%halo, div, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
    end if

    if (div_damp_order == 4) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            div2(i,j,k) = (                                                               &
              div(i+1,j,k) - 2 * div(i,j,k) + div(i-1,j,k)                                &
            ) / mesh%de_lon(j)**2 + (                                                     &
              (div(i,j+1,k) - div(i,j  ,k)) * mesh%half_cos_lat(j  ) / mesh%de_lat(j  ) - &
              (div(i,j  ,k) - div(i,j-1,k)) * mesh%half_cos_lat(j-1) / mesh%de_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block%halo, div2, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
    end if
    end associate

  end subroutine calc_div

  subroutine calc_gz_lev(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, l
    real(r8) dgz

    associate (mesh   => block%mesh      , &
               gzs    => block%static%gzs, & ! in
               t      => dstate%t        , & ! in
               ph_lev => dstate%ph_lev   , & ! in
               gz_lev => dstate%gz_lev   , & ! out
               gz     => dstate%gz       )   ! out
    do k = mesh%half_kds, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dgz = 0.0_r8
          do l = k, mesh%full_nlev
            dgz = dgz + Rd * t(i,j,l) * log(ph_lev(i,j,l+1) / ph_lev(i,j,l))
          end do
          gz_lev(i,j,k) = gzs(i,j) + dgz
        end do
      end do
    end do
    call fill_halo(block%halo, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
    ! For output
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gz(i,j,k) = 0.5_r8 * (gz_lev(i,j,k) + gz_lev(i,j,k+1))
        end do
      end do
    end do
    end associate

  end subroutine calc_gz_lev

  subroutine calc_m(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, l

    associate (mesh   => block%mesh      , &
               ph     => dstate%ph       , & ! in
               ph_lev => dstate%ph_lev   , & ! in
               gz     => dstate%gz       , & ! in
               gzs    => block%static%gzs, & ! in
               m      => dstate%m        , & ! out
               m_lon  => dstate%m_lon    , & ! out
               m_lat  => dstate%m_lat    , & ! out
               m_lev  => dstate%m_lev    , & ! out
               m_vtx  => dstate%m_vtx    )   ! out
    if (baroclinic .or. advection) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            m(i,j,k) = ph_lev(i,j,k+1) - ph_lev(i,j,k)
            if (m(i,j,k) <= 0) then
              do l = mesh%half_kds, mesh%half_kde
                print *, l, ph_lev(i,j,l)
              end do
              print *, 'phs(i,j) =', dstate%phs(i,j)
              print *, mesh%full_lon_deg(i), '(', to_str(i), ')', mesh%full_lat_deg(j), '(', to_str(j), ')', k
              print *, 'The pressure levels are not monotonic!'
              call process_stop(1)
            end if
          end do
        end do
      end do

      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            m_lev(i,j,k) = ph(i,j,k) - ph(i,j,k-1)
          end do
        end do
      end do
      ! Top boundary
      k = mesh%half_kds
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          m_lev(i,j,k) = ph(i,j,k) - ph_lev(i,j,k)
        end do
      end do
      ! Bottom boundary
      k = mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          m_lev(i,j,k) = ph_lev(i,j,k) - ph(i,j,k-1)
        end do
      end do
      call fill_halo(block%halo, m_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
    else
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          m(i,j,1) = (gz(i,j,1) - gzs(i,j)) / g
        end do
      end do
    end if

    call fill_halo(block%halo, m, full_lon=.true., full_lat=.true., full_lev=.true.)
    call average_cell_to_lon_edge(mesh, m, m_lon)
    call fill_halo(block%halo, m_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
    call average_cell_to_lat_edge(mesh, m, m_lat)
    call fill_halo(block%halo, m_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
    call interp_cell_to_vtx(mesh, m, m_vtx)
    end associate

  end subroutine calc_m

  subroutine calc_mf(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh    => block%mesh    , &
               m       => dstate%m      , & ! in
               m_lon   => dstate%m_lon  , & ! in
               m_lat   => dstate%m_lat  , & ! in
               u_lon   => dstate%u_lon  , & ! in
               v_lat   => dstate%v_lat  , & ! in
               u_lat   => dstate%u_lat  , & ! out
               v_lon   => dstate%v_lon  , & ! out
               mfx_lon => dstate%mfx_lon, & ! out
               mfy_lat => dstate%mfy_lat, & ! out
               mfy_lon => dstate%mfy_lon, & ! out
               mfx_lat => dstate%mfx_lat)   ! out
    call block%adv_batch_pt%accum_uv_cell(u_lon, v_lat, dt)
    ! call adv_calc_mass_hflx(block, block%adv_batch_pt, m, mfx_lon, mfy_lat, dt)
    ! call fill_halo(block%halo, mfx_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
    ! call fill_halo(block%halo, mfy_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%half_ids - 1, mesh%half_ide
          mfx_lon(i,j,k) = m_lon(i,j,k) * u_lon(i,j,k)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide + 1
          mfy_lat(i,j,k) = m_lat(i,j,k) * v_lat(i,j,k)
        end do
      end do
    end do
    call block%adv_batch_pt%accum_mf_cell(mfx_lon, mfy_lat)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          mfx_lat(i,j,k) = mesh%half_tangent_wgt(1,j) * (mfx_lon(i-1,j  ,k) + mfx_lon(i,j  ,k)) + &
                           mesh%half_tangent_wgt(2,j) * (mfx_lon(i-1,j+1,k) + mfx_lon(i,j+1,k))
          u_lat(i,j,k) = mfx_lat(i,j,k) / m_lat(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, u_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          mfy_lon(i,j,k) = mesh%full_tangent_wgt(1,j) * (mfy_lat(i,j-1,k) + mfy_lat(i+1,j-1,k)) + &
                           mesh%full_tangent_wgt(2,j) * (mfy_lat(i,j  ,k) + mfy_lat(i+1,j  ,k))
          v_lon(i,j,k) = mfy_lon(i,j,k) / m_lon(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, v_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine calc_mf

  subroutine calc_vor(block, dstate, u, v)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: u(block%mesh%half_ims:block%mesh%half_ime, &
                              block%mesh%full_jms:block%mesh%full_jme, &
                              block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in) :: v(block%mesh%full_ims:block%mesh%full_ime, &
                              block%mesh%half_jms:block%mesh%half_jme, &
                              block%mesh%full_kms:block%mesh%full_kme)

    integer i, j, k
    real(r8) work(dstate%mesh%half_ids:dstate%mesh%half_ide,dstate%mesh%full_nlev)
    real(r8) pole(dstate%mesh%full_nlev)

    associate (mesh  => block%mesh  , &
               u_lat => dstate%u_lat, & ! in
               vor   => dstate%vor  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%half_ids - 1, mesh%half_ide
          vor(i,j,k) = (                                                    &
            u(i  ,j,k) * mesh%de_lon(j  ) - u(i,j+1,k) * mesh%de_lon(j+1) + &
            v(i+1,j,k) * mesh%de_lat(j  ) - v(i,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            work(i,k) = -u_lat(i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            ! vor(i,j,k) = vor(i,j+1,k) / 3.0_r8 + pole(k) * 2.0_r8 / 3.0_r8
            vor(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            work(i,k) = u_lat(i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%full_nlon / mesh%area_cell(j+1)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            ! vor(i,j,k) = vor(i,j-1,k) / 3.0_r8 + pole(k) * 2.0_r8 / 3.0_r8
            vor(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    end associate

  end subroutine calc_vor

  subroutine calc_pv(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    associate (mesh  => block%mesh  , &
               m_vtx => dstate%m_vtx, & ! in
               u_lon => dstate%u_lon, & ! in
               v_lat => dstate%v_lat, & ! in
               vor   => dstate%vor  , & ! in
               pv    => dstate%pv)      ! out
    call calc_vor(block, dstate, u_lon, v_lat)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          pv(i,j,k) = (vor(i,j,k) + mesh%half_f(j)) / m_vtx(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, pv, full_lon=.false., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine calc_pv

  subroutine interp_pv_midpoint(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh   => block%mesh   , &
               pv     => dstate%pv    , & ! in
               pv_lon => dstate%pv_lon, & ! out
               pv_lat => dstate%pv_lat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          pv_lat(i,j,k) = 0.5_r8 * (pv(i-1,j,k) + pv(i,j,k))
        end do
      end do
    end do
    call fill_halo(block%halo, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          pv_lon(i,j,k) = 0.5_r8 * (pv(i,j,k) + pv(i,j-1,k))
        end do
      end do
    end do
    call fill_halo(block%halo, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)
    end associate

  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    real(r8) b
    integer i, j, k

    associate (mesh     => block%mesh   , &
               un       => dstate%u_lon , & ! in
               vn       => dstate%v_lat , & ! in
               ut       => dstate%u_lat , & ! in
               vt       => dstate%v_lon , & ! in
               pv       => dstate%pv    , & ! in
               pv_lon   => dstate%pv_lon, & ! out
               pv_lat   => dstate%pv_lat)   ! out
    select case (upwind_order_pv)
    case (1)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b  = abs(vt(i,j,k)) / (sqrt(un(i,j,k)**2 + vt(i,j,k)**2) + eps)
            pv_lon(i,j,k) = b * upwind1(sign(1.0_r8, vt(i,j,k)), upwind_wgt_pv, pv(i,j-1:j,k)) + &
                            (1 - b) * 0.5_r8 * (pv(i,j-1,k) + pv(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b  = abs(ut(i,j,k)) / (sqrt(ut(i,j,k)**2 + vn(i,j,k)**2) + eps)
            pv_lat(i,j,k) = b * upwind1(sign(1.0_r8, ut(i,j,k)), upwind_wgt_pv, pv(i-1:i,j,k)) + &
                            (1 - b) * 0.5_r8 * (pv(i-1,j,k) + pv(i,j,k))
          end do
        end do
      end do
    case (3)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b  = abs(vt(i,j,k)) / (sqrt(un(i,j,k)**2 + vt(i,j,k)**2) + eps)
            pv_lon(i,j,k) = b * upwind3(sign(1.0_r8, vt(i,j,k)), upwind_wgt_pv, pv(i,j-2:j+1,k)) + &
                            (1 - b) * 0.5_r8 * (pv(i,j-1,k) + pv(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b  = abs(ut(i,j,k)) / (sqrt(ut(i,j,k)**2 + vn(i,j,k)**2) + eps)
            pv_lat(i,j,k) = b * upwind3(sign(1.0_r8, ut(i,j,k)), upwind_wgt_pv, pv(i-2:i+1,j,k)) + &
                            (1 - b) * 0.5_r8 * (pv(i-1,j,k) + pv(i,j,k))
          end do
        end do
      end do
    end select
    call fill_halo(block%halo, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)
    call fill_halo(block%halo, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
    end associate

  end subroutine interp_pv_upwind

  subroutine interp_pv_tvd(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    real(r8) cfl
    integer i, j, k

    associate (mesh     => block%mesh   , &
               un       => dstate%u_lon , & ! in
               vn       => dstate%v_lat , & ! in
               ut       => dstate%u_lat , & ! in
               vt       => dstate%v_lon , & ! in
               pv       => dstate%pv    , & ! in
               pv_lon   => dstate%pv_lon, & ! out
               pv_lat   => dstate%pv_lat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          cfl = vt(i,j,k) * dt / mesh%le_lon(j)
          pv_lon(i,j,k) = tvd(cfl, pv(i,j-2,k), pv(i,j-1,k), pv(i,j,k), pv(i,j+1,k))
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          cfl = ut(i,j,k) * dt / mesh%le_lat(j)
          pv_lat(i,j,k) = tvd(cfl, pv(i-2,j,k), pv(i-1,j,k), pv(i,j,k), pv(i+1,j,k))
        end do
      end do
    end do
    call fill_halo(block%halo, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)
    call fill_halo(block%halo, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
    end associate

  end subroutine interp_pv_tvd

  subroutine calc_coriolis(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh    => block%mesh    , &
               mfx_lon => dstate%mfx_lon, & ! in
               mfy_lat => dstate%mfy_lat, & ! in
               mfy_lon => dstate%mfy_lon, & ! in
               mfx_lat => dstate%mfx_lat, & ! in
               pv_lon  => dstate%pv_lon , & ! in
               pv_lat  => dstate%pv_lat , & ! in
               qhu     => dtend%qhu     , & ! out
               qhv     => dtend%qhv     )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          if (coriolis_scheme == 1) then
            qhu(i,j,k) = (                                                 &
              mesh%half_tangent_wgt(1,j) * (                               &
                mfx_lon(i-1,j  ,k) * (pv_lat(i,j,k) + pv_lon(i-1,j  ,k)) + &
                mfx_lon(i  ,j  ,k) * (pv_lat(i,j,k) + pv_lon(i  ,j  ,k))   &
              ) +                                                          &
              mesh%half_tangent_wgt(2,j) * (                               &
                mfx_lon(i-1,j+1,k) * (pv_lat(i,j,k) + pv_lon(i-1,j+1,k)) + &
                mfx_lon(i  ,j+1,k) * (pv_lat(i,j,k) + pv_lon(i  ,j+1,k))   &
              )                                                            &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            qhu(i,j,k) = mfx_lat(i,j,k) * pv_lat(i,j,k)
          end if
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          if (coriolis_scheme == 1) then
            qhv(i,j,k) = (                                                 &
              mesh%full_tangent_wgt(1,j) * (                               &
                mfy_lat(i  ,j-1,k) * (pv_lon(i,j,k) + pv_lat(i  ,j-1,k)) + &
                mfy_lat(i+1,j-1,k) * (pv_lon(i,j,k) + pv_lat(i+1,j-1,k))   &
              ) +                                                          &
              mesh%full_tangent_wgt(2,j) * (                               &
                mfy_lat(i  ,j  ,k) * (pv_lon(i,j,k) + pv_lat(i  ,j  ,k)) + &
                mfy_lat(i+1,j  ,k) * (pv_lon(i,j,k) + pv_lat(i+1,j  ,k))   &
              )                                                            &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            qhv(i,j,k) = mfy_lon(i,j,k) * pv_lon(i,j,k)
          end if
        end do
      end do
    end do
    end associate

  end subroutine calc_coriolis

  subroutine calc_grad_ke(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh    => block%mesh   , &
               ke      => dstate%ke    , & ! in
               dkedlon => dtend%dkedlon, & ! out
               dkedlat => dtend%dkedlat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          dkedlon(i,j,k) = (ke(i+1,j,k) - ke(i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          dkedlat(i,j,k) = (ke(i,j+1,k) - ke(i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine calc_grad_ke

  subroutine calc_grad_mf(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) work(dstate%mesh%full_ids:dstate%mesh%full_ide,dstate%mesh%full_nlev)
    real(r8) pole(dstate%mesh%full_nlev)

    associate (mesh    => block%mesh    , &
               mfx_lon => dstate%mfx_lon, & ! in
               mfy_lat => dstate%mfy_lat, & ! in
               dmfdlon => dtend%dmfdlon , & ! out
               dmfdlat => dtend%dmfdlat )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dmfdlon(i,j,k) = (                  &
            mfx_lon(i,j,k) - mfx_lon(i-1,j,k) &
          ) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dmfdlat(i,j,k) = (                      &
            mfy_lat(i,j  ,k) * mesh%le_lat(j  ) - &
            mfy_lat(i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy_lat(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dmfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -mfy_lat(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dmfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
    end associate

  end subroutine calc_grad_mf

  subroutine calc_grad_ptf(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) work(dstate%mesh%full_ids:dstate%mesh%full_ide,dstate%mesh%full_nlev)
    real(r8) pole(dstate%mesh%full_nlev)

    associate (mesh     => block%mesh    , &
               gz_lev   => dstate%gz_lev , & ! in
               ph       => dstate%ph     , & ! in
               pt       => dstate%pt     , & ! in
               ptf_lon  => dstate%ptf_lon, & ! out
               ptf_lat  => dstate%ptf_lat, & ! out
               ptf_lev  => dstate%ptf_lev, & ! out
               dptfdlon => dtend%dptfdlon, & ! out
               dptfdlat => dtend%dptfdlat, & ! out
               dptfdlev => dtend%dptfdlev)   ! out
    call adv_calc_tracer_hflx(block, block%adv_batch_pt, pt, ptf_lon, ptf_lat, dt)
    call fill_halo(block%halo, ptf_lon, full_lon=.false., full_lat=.true., full_lev=.true., &
                   south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(block%halo, ptf_lat, full_lon=.true., full_lat=.false., full_lev=.true., &
                   north_halo=.false.,  west_halo=.false., east_halo=.false.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dptfdlon(i,j,k) = (                 &
            ptf_lon(i,j,k) - ptf_lon(i-1,j,k) &
          ) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dptfdlat(i,j,k) = (                     &
            ptf_lat(i,j  ,k) * mesh%le_lat(j  ) - &
            ptf_lat(i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = ptf_lat(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dptfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -ptf_lat(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dptfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
    ! --------------------------------- FFSL -----------------------------------
    ! Set upper and lower boundary conditions.
    do k = mesh%full_kds - 1, mesh%full_kds - 2, -1
      pt(:,:,k) = 2 * pt(:,:,k+1) - pt(:,:,k+2)
    end do
    do k = mesh%full_kde + 1, mesh%full_kde + 2
      pt(:,:,k) = 2 * pt(:,:,k-1) - pt(:,:,k-2)
    end do
    call adv_calc_tracer_vflx(block, block%adv_batch_pt, pt, ptf_lev, dt)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dptfdlev(i,j,k) = ptf_lev(i,j,k+1) - ptf_lev(i,j,k)
        end do
      end do
    end do
    ! ------------------------------ Taylor 2020 -------------------------------
    ! do k = mesh%half_kds + 1, mesh%half_kde - 1
    !   do j = mesh%full_jds, mesh%full_jde
    !     do i = mesh%full_ids, mesh%full_ide
    !       ptf_lev(i,j,k) = -(gz_lev(i,j,k+1) - gz_lev(i,j,k-1)) / (2 * cpd) / &
    !         (ph(i,j,k)**Rd_o_cpd - ph(i,j,k-1)**Rd_o_cpd) * p0**Rd_o_cpd
    !     end do
    !   end do
    ! end do
    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds, mesh%full_jde
    !     do i = mesh%full_ids, mesh%full_ide
    !       dptfdlev(i,j,k) = dstate%we_lev(i,j,k+1) * ptf_lev(i,j,k+1) - &
    !                         dstate%we_lev(i,j,k  ) * ptf_lev(i,j,k  )
    !     end do
    !   end do
    ! end do
    end associate

  end subroutine calc_grad_ptf

  subroutine calc_dphsdt(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh    => block%mesh   , &
               dmfdlon => dtend%dmfdlon, & ! in
               dmfdlat => dtend%dmfdlat, & ! in
               dphs    => dtend%dphs   )   ! out
    dtend%dphs = 0
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dphs(i,j) = dphs(i,j) - dmfdlon(i,j,k) - dmfdlat(i,j,k)
        end do
      end do
    end do
    end associate

  end subroutine calc_dphsdt

  subroutine calc_wedudlev_wedvdlev(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    ! Follow SB81 vertical advection discretization.

    associate (mesh       => block%mesh       , &
               u          => dstate%u_lon     , & ! in
               v          => dstate%v_lat     , & ! in
               m_lon      => dstate%m_lon     , & ! in
               m_lat      => dstate%m_lat     , & ! in
               we_lev_lon => dstate%we_lev_lon, & ! in
               we_lev_lat => dstate%we_lev_lat, & ! in
               wedudlev   => dtend%wedudlev   , & ! out
               wedvdlev   => dtend%wedvdlev   )   ! out
    do k = mesh%full_kds + 1, mesh%full_kde - 1
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          wedudlev(i,j,k) = (                               &
            we_lev_lon(i,j,k+1) * (u(i,j,k+1) - u(i,j,k)) + &
            we_lev_lon(i,j,k  ) * (u(i,j,k) - u(i,j,k-1))   &
          ) / m_lon(i,j,k) / 2.0_r8
        end do
      end do
    end do
    k = mesh%full_kds
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids, mesh%half_ide
        wedudlev(i,j,k) = (we_lev_lon(i,j,k+1) * &
          (u(i,j,k+1) - u(i,j,k))) / m_lon(i,j,k) / 2.0_r8
      end do
    end do
    k = mesh%full_kde
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%half_ids, mesh%half_ide
        wedudlev(i,j,k) = (we_lev_lon(i,j,k  ) * &
          (u(i,j,k) - u(i,j,k-1))) / m_lon(i,j,k) / 2.0_r8
      end do
    end do

    do k = mesh%full_kds + 1, mesh%full_kde - 1
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          wedvdlev(i,j,k) = (                                 &
            we_lev_lat(i,j,k+1) * (v(i,j,k+1) - v(i,j,k  )) + &
            we_lev_lat(i,j,k  ) * (v(i,j,k  ) - v(i,j,k-1))   &
          ) / m_lat(i,j,k) / 2.0_r8
        end do
      end do
    end do
    k = mesh%full_kds
    do j = mesh%half_jds, mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        wedvdlev(i,j,k) = (we_lev_lat(i,j,k+1) * &
          (v(i,j,k+1) - v(i,j,k))) / m_lat(i,j,k) / 2.0_r8
      end do
    end do
    k = mesh%full_kde
    do j = mesh%half_jds, mesh%half_jde
      do i = mesh%full_ids, mesh%full_ide
        wedvdlev(i,j,k) = (we_lev_lat(i,j,k  ) * &
          (v(i,j,k) - v(i,j,k-1))) / m_lat(i,j,k) / 2.0_r8
      end do
    end do
    end associate

  end subroutine calc_wedudlev_wedvdlev

end module operators_mod
