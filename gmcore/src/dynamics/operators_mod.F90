module operators_mod

  use const_mod
  use perf_mod
  use vert_coord_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod, only: proc, process_stop
  use formula_mod
  use namelist_mod
  use tracer_mod
  use log_mod
  use pgf_mod
  use adv_mod
  use interp_mod
  use filter_mod

  implicit none

  private

  public operators_init
  public operators_prepare
  public calc_mg
  public calc_ph
  public calc_p
  public calc_omg
  public calc_dmg
  public calc_t
  public calc_rhod
  public calc_gz_lev
  public calc_we_lev
  public calc_div
  public calc_vor
  public calc_coriolis
  public calc_grad_ke
  public calc_grad_mf
  public calc_grad_ptf
  public calc_dmgsdt
  public calc_wedudlev_wedvdlev

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

  interface
    subroutine interp_pv_interface(block, dstate, dt, substep)
      import block_type, dstate_type, r8
      type(block_type), intent(inout) :: block
      type(dstate_type), intent(inout) :: dstate
      real(r8), intent(in) :: dt
      integer, intent(in) :: substep
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
      if (baroclinic    ) call calc_mg    (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_dmg                       (blocks(iblk), blocks(iblk)%dstate(itime))
      if (baroclinic    ) call calc_ph    (blocks(iblk), blocks(iblk)%dstate(itime))
      if (baroclinic    ) call calc_t     (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_mf                        (blocks(iblk), blocks(iblk)%dstate(itime), dt)
      call calc_ke                        (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_pv                        (blocks(iblk), blocks(iblk)%dstate(itime))
      call interp_pv                      (blocks(iblk), blocks(iblk)%dstate(itime), dt, total_substeps)
      if (baroclinic    ) call calc_gz_lev(blocks(iblk), blocks(iblk)%dstate(itime))
      if (baroclinic    ) call calc_rhod  (blocks(iblk), blocks(iblk)%dstate(itime))
      if (nonhydrostatic) call calc_p     (blocks(iblk), blocks(iblk)%dstate(itime))
      call pgf_prepare                    (blocks(iblk), blocks(iblk)%dstate(itime))
      call tracer_calc_qm                 (blocks(iblk))
      if (nonhydrostatic) call fill_halo(blocks(iblk)%dstate(itime)%gz_lev)
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, dstate, dt, pass, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: substep

    select case (pass)
    ! --------------------------------------------------------------------------
    case (all_pass)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate)
      call calc_pv                        (block, dstate)
      call interp_pv                      (block, dstate, dt, substep)
      if (baroclinic    ) call calc_t     (block, dstate)
      if (hydrostatic   ) call calc_gz_lev(block, dstate)
      if (hydrostatic   ) call calc_rhod  (block, dstate)
      call pgf_prepare                    (block, dstate)
    ! --------------------------------------------------------------------------
    case (forward_pass)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate)
      call calc_pv                        (block, dstate)
      call interp_pv                      (block, dstate, dt, substep)
    ! --------------------------------------------------------------------------
    case (backward_pass)
      if (baroclinic    ) call calc_t     (block, dstate)
      if (hydrostatic   ) call calc_gz_lev(block, dstate)
      if (hydrostatic   ) call calc_rhod  (block, dstate)
      call pgf_prepare                    (block, dstate)
    end select

  end subroutine operators_prepare_2

  subroutine calc_mg(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_mg')

    associate (mesh    => block%mesh    , &
               mgs     => dstate%mgs    , & ! in
               mg_lev  => dstate%mg_lev , & ! out
               mg      => dstate%mg     )   ! out
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          mg_lev%d(i,j,k) = vert_coord_calc_mg_lev(k, mgs%d(i,j), block%static%ref_ps_perb%d(i,j))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          mg%d(i,j,k) = 0.5_r8 * (mg_lev%d(i,j,k) + mg_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_mg')

  end subroutine calc_mg

  subroutine calc_ph(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_ph')

    associate (mesh    => block%mesh          , &
               mg_lev  => dstate%mg_lev       , & ! in
               dmg     => dstate%dmg          , & ! in
               qm      => tracers(block%id)%qm, & ! in
               ph_lev  => dstate%ph_lev       , & ! out
               pkh_lev => block%aux%pkh_lev   , & ! out
               ph      => dstate%ph           , & ! out
               phs     => dstate%phs          , & ! pointer
               ps      => dstate%ps           )   ! out
    k = mesh%half_kds
    ph_lev%d(:,:,k) = mg_lev%d(:,:,k)
    pkh_lev%d(:,:,k) = ph_lev%d(:,:,k)**rd_o_cpd
    do k = mesh%half_kds + 1, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          ph_lev%d(i,j,k) = ph_lev%d(i,j,k-1) + dmg%d(i,j,k-1) * (1 + qm%d(i,j,k-1))
          pkh_lev%d(i,j,k) = ph_lev%d(i,j,k)**rd_o_cpd
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          ph%d(i,j,k) = 0.5_r8 * (ph_lev%d(i,j,k) + ph_lev%d(i,j,k+1))
        end do
      end do
    end do
    ! NOTE: Move this to other place?
    if (hydrostatic) ps%d = phs%d
    end associate

    call perf_stop('calc_ph')

  end subroutine calc_ph

  subroutine calc_p(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    real(r8), parameter :: p0 = 1.0e5_r8
    integer i, j, k

    call perf_start('calc_p')

    associate (mesh  => block%mesh  , & ! in
               rhod  => dstate%rhod , & ! in
               pt    => dstate%pt   , & ! in
               p     => dstate%p    , & ! out
               p_lev => dstate%p_lev)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          p%d(i,j,k) = p0 * (Rd * pt%d(i,j,k) * rhod%d(i,j,k) / p0)**cpd_o_cvd
        end do
      end do
    end do
    call interp_run(p, p_lev)
    p_lev%d(:,:,1) = ptop
    call fill_halo(p_lev)
    end associate

    call perf_stop('calc_p')

  end subroutine calc_p

  subroutine calc_omg(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) sum_dmf(block%mesh%full_ids:block%mesh%full_ide, &
                     block%mesh%full_jds:block%mesh%full_jde)

    associate (mesh  => block%mesh   , &
               ph    => dstate%ph    , &
               u_lon => dstate%u_lon , &
               v_lat => dstate%v_lat , &
               dmf   => block%aux%dmf, &
               div   => block%aux%div, &
               omg   => block%aux%omg)
    sum_dmf = 0
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k)
          omg%d(i,j,k) = 0.5_r8 * ((                                              &
            u_lon%d(i  ,j,k) * (ph%d(i,j,k) + ph%d(i+1,j,k)) -                    &
            u_lon%d(i-1,j,k) * (ph%d(i,j,k) + ph%d(i-1,j,k))                      &
          ) * mesh%le_lon(j) - (                                                  &
            v_lat%d(i,j  ,k) * (ph%d(i,j,k) + ph%d(i,j+1,k)) * mesh%le_lat(j  ) - &
            v_lat%d(i,j-1,k) * (ph%d(i,j,k) + ph%d(i,j-1,k)) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j) - ph%d(i,j,k) * div%d(i,j,k) - sum_dmf(i,j)
        end do
      end do
    end do
    end associate

  end subroutine calc_omg

  subroutine calc_t(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_t')

    associate (mesh => block%mesh         , &
               pt   => dstate%pt          , & ! in
               ph   => dstate%ph          , & ! in
               q    => tracers(block%id)%q, & ! in
               t    => dstate%t           , & ! out
               tv   => dstate%tv          )   ! out
    if (idx_qv > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            t%d(i,j,k) = temperature(pt%d(i,j,k), ph%d(i,j,k), q%d(i,j,k,idx_qv))
            tv%d(i,j,k) = virtual_temperature_from_modified_potential_temperature(pt%d(i,j,k), ph%d(i,j,k)**rd_o_cpd, q%d(i,j,k,idx_qv))
          end do
        end do
      end do
    else
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            t%d(i,j,k) = temperature(pt%d(i,j,k), ph%d(i,j,k), 0.0_r8)
            tv%d(i,j,k) = t%d(i,j,k)
          end do
        end do
      end do
    end if
    end associate

    call perf_stop('calc_t')

  end subroutine calc_t

  subroutine calc_rhod(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_rhod')

    associate (mesh   => block%mesh   , &
               gz_lev => dstate%gz_lev, & ! in
               dmg    => dstate%dmg   , & ! in
               rhod   => dstate%rhod  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          rhod%d(i,j,k) = dmg%d(i,j,k) / (gz_lev%d(i,j,k) - gz_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_rhod')

  end subroutine calc_rhod

  subroutine calc_we_lev(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(in) :: dtend
    real(r8), intent(in) :: dt
    
    integer i, j, k
    integer result
    real(r8) sum_dmf(block%mesh%full_ids:block%mesh%full_ide, &
                     block%mesh%full_jds:block%mesh%full_jde)


    call perf_start('calc_we_lev')

    associate (mesh       => block%mesh          , &
               dmf        => block%aux%dmf       , & ! in
               dmgs       => dtend%dmgs          , & ! in
               we_lev     => dstate%we_lev       , & ! out
               we_lev_lon => block%aux%we_lev_lon, & ! out
               we_lev_lat => block%aux%we_lev_lat)   ! out


    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          if (k .eq. mesh%half_kds + 1) then 
            sum_dmf(i,j) = sum(dmf%d(i,j,1:k-1))
          else 
            sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k-1)
          end if
          ! we_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum(dmf%d(i,j,1:k-1))
          we_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum_dmf(i,j)
        end do
      end do
    end do
    call fill_halo(we_lev, west_halo=.false., south_halo=.false.)

    call interp_run(we_lev, we_lev_lon)
    call interp_run(we_lev, we_lev_lat)
    end associate

    call perf_stop('calc_we_lev')

  end subroutine calc_we_lev

  subroutine calc_ke(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) ke_vtx(4)
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call perf_start('calc_ke')

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, & ! in
               v    => dstate%v_lat, & ! in
               ke   => block%aux%ke)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          ke%d(i,j,k) = (mesh%area_lon_west (j  ) * u%d(i-1,j  ,k)**2 + &
                         mesh%area_lon_east (j  ) * u%d(i  ,j  ,k)**2 + &
                         mesh%area_lat_north(j-1) * v%d(i  ,j-1,k)**2 + &
                         mesh%area_lat_south(j  ) * v%d(i  ,j  ,k)**2   &
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
            ke_vtx(1) = (                                    &
              mesh%area_lat_east (j  ) * v%d(i-1,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v%d(i  ,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u%d(i-1,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u%d(i-1,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke_vtx(2) = (                                    &
              mesh%area_lat_east (j-1) * v%d(i-1,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v%d(i  ,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u%d(i-1,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u%d(i-1,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(3) = (                                    &
              mesh%area_lat_east (j-1) * v%d(i  ,j-1,k)**2 + &
              mesh%area_lat_west (j-1) * v%d(i+1,j-1,k)**2 + &
              mesh%area_lon_north(j-1) * u%d(i  ,j-1,k)**2 + &
              mesh%area_lon_south(j  ) * u%d(i  ,j  ,k)**2   &
            ) / mesh%area_vtx(j-1)
            ke_vtx(4) = (                                    &
              mesh%area_lat_east (j  ) * v%d(i  ,j  ,k)**2 + &
              mesh%area_lat_west (j  ) * v%d(i+1,j  ,k)**2 + &
              mesh%area_lon_north(j  ) * u%d(i  ,j  ,k)**2 + &
              mesh%area_lon_south(j+1) * u%d(i  ,j+1,k)**2   &
            ) / mesh%area_vtx(j)
            ke%d(i,j,k) = (1.0_r8 - ke_cell_wgt) * (             &
              (ke_vtx(1) + ke_vtx(4)) * mesh%area_subcell(2,j) + &
              (ke_vtx(2) + ke_vtx(3)) * mesh%area_subcell(1,j)   &
            ) / mesh%area_cell(j) + ke_cell_wgt * ke%d(i,j,k)
          end do
        end do
      end do
    end if

    ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = 0.5_r8 * (v%d(i,j,k)**2 + block%aux%u_lat%d(i,j,k)**2)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = 0.5_r8 * (v%d(i,j-1,k)**2 + block%aux%u_lat%d(i,j-1,k)**2)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke%d(i,j,k) = pole(k)
        end do
      end do
    end if
    end associate

    call perf_stop('calc_ke')

  end subroutine calc_ke

  subroutine calc_div(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)
    integer i, j, k

    call perf_start('calc_div')

    associate (mesh => block%mesh    , &
               u    => dstate%u_lon  , & ! in
               v    => dstate%v_lat  , & ! in
               div  => block%aux%div , & ! out
               div2 => block%aux%div2)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          div%d(i,j,k) = (                                                      &
            (u%d(i,j,k) * mesh%le_lon(j) - u%d(i-1,  j,k) * mesh%le_lon(j  )) + &
            (v%d(i,j,k) * mesh%le_lat(j) - v%d(i  ,j-1,k) * mesh%le_lat(j-1))   &
          ) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -v%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          div%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (div_damp_order == 4) then
      call fill_halo(div)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            div2%d(i,j,k) = (                                                                 &
              div%d(i+1,j,k) - 2 * div%d(i,j,k) + div%d(i-1,j,k)                              &
            ) / mesh%de_lon(j)**2 + (                                                         &
              (div%d(i,j+1,k) - div%d(i,j  ,k)) * mesh%half_cos_lat(j  ) / mesh%de_lat(j  ) - &
              (div%d(i,j  ,k) - div%d(i,j-1,k)) * mesh%half_cos_lat(j-1) / mesh%de_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
          end do
        end do
      end do
    end if
    end associate

    call perf_stop('calc_div')

  end subroutine calc_div

  subroutine calc_gz_lev(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_gz_lev')

    associate (mesh   => block%mesh      , &
               gzs    => block%static%gzs, & ! in
               tv     => dstate%tv       , & ! in
               ph_lev => dstate%ph_lev   , & ! in
               gz_lev => dstate%gz_lev   , & ! out
               gz     => dstate%gz       )   ! out
    do k = mesh%half_kde - 1, mesh%half_kds, -1
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          gz_lev%d(i,j,k) = gz_lev%d(i,j,k+1) + rd * tv%d(i,j,k) * log(ph_lev%d(i,j,k+1) / ph_lev%d(i,j,k))
        end do
      end do
    end do
    ! For output
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          gz%d(i,j,k) = 0.5_r8 * (gz_lev%d(i,j,k) + gz_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_gz_lev')

  end subroutine calc_gz_lev

  subroutine calc_dmg(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, l

    call perf_start('calc_dmg')

    associate (mesh    => block%mesh       , &
               mg      => dstate%mg        , & ! in
               mg_lev  => dstate%mg_lev    , & ! in
               gz      => dstate%gz        , & ! in
               gzs     => block%static%gzs , & ! in
               dmg     => dstate%dmg       , & ! out
               dmg_lon => block%aux%dmg_lon, & ! out
               dmg_lat => block%aux%dmg_lat, & ! out
               dmg_lev => dstate%dmg_lev   , & ! out
               dmg_vtx => block%aux%dmg_vtx)   ! out
    if (baroclinic .or. advection) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dmg%d(i,j,k) = mg_lev%d(i,j,k+1) - mg_lev%d(i,j,k)
            if (dmg%d(i,j,k) <= 0) then
              do l = mesh%half_kds, mesh%half_kde
                print *, l, mg_lev%d(i,j,l)
              end do
              print *, 'mgs(i,j) =', dstate%mgs%d(i,j)
              print *, mesh%full_lon_deg(i), '(', to_str(i), ')', mesh%full_lat_deg(j), '(', to_str(j), ')', k
              call log_warning('The dry-air weight levels are not monotonic!', __FILE__, __LINE__)
              call process_stop(1)
            end if
          end do
        end do
      end do

      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dmg_lev%d(i,j,k) = mg%d(i,j,k) - mg%d(i,j,k-1)
          end do
        end do
      end do
      ! Top boundary
      k = mesh%half_kds
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dmg_lev%d(i,j,k) = mg%d(i,j,k) - mg_lev%d(i,j,k)
        end do
      end do
      ! Bottom boundary
      k = mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dmg_lev%d(i,j,k) = mg_lev%d(i,j,k) - mg%d(i,j,k-1)
        end do
      end do
      call fill_halo(dmg_lev)
    else
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dmg%d(i,j,1) = (gz%d(i,j,1) - gzs%d(i,j)) / g
        end do
      end do
    end if

    call fill_halo(dmg)
    call average_run(dmg, dmg_lon)
    call fill_halo(dmg_lon)
    call average_run(dmg, dmg_lat)
    call fill_halo(dmg_lat)
    call interp_run(dmg, dmg_vtx)
    end associate

    call perf_stop('calc_dmg')

  end subroutine calc_dmg

  subroutine calc_mf(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    call perf_start('calc_mf')

    associate (mesh    => block%mesh       , &
               dmg     => dstate%dmg       , & ! in
               dmg_lon => block%aux%dmg_lon, & ! in
               dmg_lat => block%aux%dmg_lat, & ! in
               u_lon   => dstate%u_lon     , & ! in
               v_lat   => dstate%v_lat     , & ! in
               u_lat   => block%aux%u_lat  , & ! out
               v_lon   => block%aux%v_lon  , & ! out
               mfx_lon => block%aux%mfx_lon, & ! out
               mfy_lat => block%aux%mfy_lat, & ! out
               mfy_lon => block%aux%mfy_lon, & ! out
               mfx_lat => block%aux%mfx_lat)   ! out
    
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%half_ids - 1, mesh%half_ide
          mfx_lon%d(i,j,k) = dmg_lon%d(i,j,k) * u_lon%d(i,j,k)
        end do
      end do


    ! do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          mfx_lat%d(i,j,k) = block%static%tg_wgt_lat(1,j) * (mfx_lon%d(i-1,j  ,k) + mfx_lon%d(i,j  ,k)) + &
                             block%static%tg_wgt_lat(2,j) * (mfx_lon%d(i-1,j+1,k) + mfx_lon%d(i,j+1,k))
          u_lat%d(i,j,k) = mfx_lat%d(i,j,k) / dmg_lat%d(i,j,k)
        end do
      end do

      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide + 1
          mfy_lat%d(i,j,k) = dmg_lat%d(i,j,k) * v_lat%d(i,j,k)
        end do
      end do
    ! end do
    ! do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          mfy_lon%d(i,j,k) = block%static%tg_wgt_lon(1,j) * (mfy_lat%d(i,j-1,k) + mfy_lat%d(i+1,j-1,k)) + &
                             block%static%tg_wgt_lon(2,j) * (mfy_lat%d(i,j  ,k) + mfy_lat%d(i+1,j  ,k))
          v_lon%d(i,j,k) = mfy_lon%d(i,j,k) / dmg_lon%d(i,j,k)
        end do
      end do
    end do


    call fill_halo(u_lat)

    ! do k = mesh%full_kds, mesh%full_kde

    ! end do
    call fill_halo(v_lon)
    end associate

    call perf_stop('calc_mf')

  end subroutine calc_mf

  subroutine calc_vor(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) work(block%mesh%half_ids:block%mesh%half_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call perf_start('calc_vor')

    associate (mesh  => block%mesh     , &
               u_lon => dstate%u_lon   , & ! in
               v_lat => dstate%v_lat   , & ! in
               u_lat => block%aux%u_lat, & ! in
               vor   => block%aux%vor  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          vor%d(i,j,k) = (                                                            &
            u_lon%d(i  ,j,k) * mesh%de_lon(j) - u_lon%d(i,j+1,k) * mesh%de_lon(j+1) + &
            v_lat%d(i+1,j,k) * mesh%de_lat(j) - v_lat%d(i,j  ,k) * mesh%de_lat(j  )   &
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
            work(i,k) = -u_lat%d(i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%full_nlon / mesh%area_cell(j)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            vor%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            work(i,k) = u_lat%d(i,j,k) * mesh%le_lat(j)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%full_nlon / mesh%area_cell(j+1)
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            vor%d(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    end associate

    call perf_stop('calc_vor')

  end subroutine calc_vor

  subroutine calc_pv(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_pv')

    associate (mesh    => block%mesh       , &
               dmg_vtx => block%aux%dmg_vtx, & ! in
               vor     => block%aux%vor    , & ! in
               pv      => block%aux%pv     )   ! out
    call calc_vor(block, dstate)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          pv%d(i,j,k) = (vor%d(i,j,k) + block%static%f_lat(j)) / dmg_vtx%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pv)
    end associate

    call perf_stop('calc_pv')

  end subroutine calc_pv

  subroutine interp_pv_midpoint(block, dstate, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    integer i, j, k

    associate (mesh   => block%mesh      , &
               pv     => block%aux%pv    , & ! in
               pv_lon => block%aux%pv_lon, & ! out
               pv_lat => block%aux%pv_lat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          pv_lat%d(i,j,k) = 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          pv_lon%d(i,j,k) = 0.5_r8 * (pv%d(i,j,k) + pv%d(i,j-1,k))
        end do
      end do
    end do
    call fill_halo(pv_lon, east_halo=.false., south_halo=.false.)
    call fill_halo(pv_lat, west_halo=.false., north_halo=.false.)
    end associate

  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, dstate, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    real(r8) b
    integer i, j, k

    ! if (substep < total_substeps) then
    !   call interp_pv_midpoint(block, dstate, dt, substep)
    !   return
    ! end if

    call perf_start('interp_pv_upwind')

    associate (mesh   => block%mesh      , &
               un     => dstate%u_lon    , & ! in
               vn     => dstate%v_lat    , & ! in
               ut     => block%aux%u_lat , & ! in
               vt     => block%aux%v_lon , & ! in
               pv     => block%aux%pv    , & ! in
               pv_lon => block%aux%pv_lon, & ! out
               pv_lat => block%aux%pv_lat)   ! out
    select case (upwind_order_pv)
    case (1)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind1(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-1:j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind1(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-1:i,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (3)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind3(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-2:j+1,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b  = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind3(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-2:i+1,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (5)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind5(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-3:j+2,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind5(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-3:i+2,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    end select
    call fill_halo(pv_lon, east_halo=.false., south_halo=.false.)
    call fill_halo(pv_lat, west_halo=.false., north_halo=.false.)
    end associate

    call perf_stop('interp_pv_upwind')

  end subroutine interp_pv_upwind

  subroutine calc_coriolis(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k

    call perf_start('calc_coriolis')

    associate (mesh    => block%mesh       , &
               mfx_lon => block%aux%mfx_lon, & ! in
               mfy_lat => block%aux%mfy_lat, & ! in
               mfy_lon => block%aux%mfy_lon, & ! in
               mfx_lat => block%aux%mfx_lat, & ! in
               pv_lon  => block%aux%pv_lon , & ! in
               pv_lat  => block%aux%pv_lat , & ! in
               du      => dtend%du         , & ! out
               dv      => dtend%dv         )   ! out
    select case (coriolis_scheme)
    case (1)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            tmp = - (                                                            &
              block%static%tg_wgt_lat(1,j) * (                                   &
                mfx_lon%d(i-1,j  ,k) * (pv_lat%d(i,j,k) + pv_lon%d(i-1,j  ,k)) + &
                mfx_lon%d(i  ,j  ,k) * (pv_lat%d(i,j,k) + pv_lon%d(i  ,j  ,k))   &
              ) +                                                                &
              block%static%tg_wgt_lat(2,j) * (                                   &
                mfx_lon%d(i-1,j+1,k) * (pv_lat%d(i,j,k) + pv_lon%d(i-1,j+1,k)) + &
                mfx_lon%d(i  ,j+1,k) * (pv_lat%d(i,j,k) + pv_lon%d(i  ,j+1,k))   &
              )                                                                  &
            ) * 0.5_r8
            dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dvdt_coriolis%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            tmp = (                                                              &
              block%static%tg_wgt_lon(1,j) * (                                   &
                mfy_lat%d(i  ,j-1,k) * (pv_lon%d(i,j,k) + pv_lat%d(i  ,j-1,k)) + &
                mfy_lat%d(i+1,j-1,k) * (pv_lon%d(i,j,k) + pv_lat%d(i+1,j-1,k))   &
              ) +                                                                &
              block%static%tg_wgt_lon(2,j) * (                                   &
                mfy_lat%d(i  ,j  ,k) * (pv_lon%d(i,j,k) + pv_lat%d(i  ,j  ,k)) + &
                mfy_lat%d(i+1,j  ,k) * (pv_lon%d(i,j,k) + pv_lat%d(i+1,j  ,k))   &
              )                                                                  &
            ) * 0.5_r8
            du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
            dtend%dudt_coriolis%d(i,j,k) = tmp
#endif
          end do
        end do
      end do
    case (2)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            dv%d(i,j,k) = dv%d(i,j,k) - mfx_lat%d(i,j,k) * pv_lat%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            du%d(i,j,k) = du%d(i,j,k) + mfy_lon%d(i,j,k) * pv_lon%d(i,j,k)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('calc_coriolis')

  end subroutine calc_coriolis

  subroutine calc_grad_ke(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k

    call perf_start('calc_grad_ke')

    associate (mesh => block%mesh  , &
               ke   => block%aux%ke, & ! in
               du   => dtend%du    , & ! out
               dv   => dtend%dv    )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp = -(ke%d(i+1,j,k) - ke%d(i,j,k)) / mesh%de_lon(j)
          du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_dkedx%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp = -(ke%d(i,j+1,k) - ke%d(i,j,k)) / mesh%de_lat(j)
          dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_dkedy%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_grad_ke')

  end subroutine calc_grad_ke

  subroutine calc_grad_mf(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call perf_start('calc_grad_mf')

    associate (mesh    => block%mesh       , &
               mfx_lon => block%aux%mfx_lon, & ! in
               mfy_lat => block%aux%mfy_lat, & ! in
               dmf     => block%aux%dmf    )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dmf%d(i,j,k) = ((                         &
            mfx_lon%d(i,j,k) - mfx_lon%d(i-1,j,k)   &
          ) * mesh%le_lon(j) + (                    &
            mfy_lat%d(i,j  ,k) * mesh%le_lat(j  ) - &
            mfy_lat%d(i,j-1,k) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy_lat%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dmf%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = -mfy_lat%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dmf%d(i,j,k) = pole(k)
        end do
      end do
    end if
    end associate

    call perf_stop('calc_grad_mf')

  end subroutine calc_grad_mf

  subroutine calc_grad_ptf(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call perf_start('calc_grad_ptf')

    associate (mesh     => block%filter_mesh, &
               u_lon    => dstate%u_lon     , & ! in
               v_lat    => dstate%v_lat     , & ! in
               we_lev   => dstate%we_lev    , & ! in
               mfx_lon  => block%aux%mfx_lon, & ! in
               mfy_lat  => block%aux%mfy_lat, & ! in
               dmg_lev  => dstate%dmg_lev   , & ! in
               pt       => dstate%pt        , & ! in
               ptf_lon  => block%aux%ptf_lon, & ! out
               ptf_lat  => block%aux%ptf_lat, & ! out
               ptf_lev  => block%aux%ptf_lev, & ! out
               dpt      => dtend%dpt        )   ! out
    call block%adv_batch_pt%set_wind(u_lon, v_lat, we_lev, mfx_lon, mfy_lat, dmg_lev)
    call adv_calc_tracer_hflx(block%adv_batch_pt, pt, ptf_lon, ptf_lat, dt)
    call fill_halo(ptf_lon, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(ptf_lat, north_halo=.false.,  west_halo=.false., east_halo=.false.)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          dpt%d(i,j,k) = -((                        &
            ptf_lon%d(i,j,k) - ptf_lon%d(i-1,j,k)   &
          ) * mesh%le_lon(j) + (                    &
            ptf_lat%d(i,j  ,k) * mesh%le_lat(j  ) - &
            ptf_lat%d(i,j-1,k) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = ptf_lat%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dpt%d(i,j,k) = -pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = ptf_lat%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          dpt%d(i,j,k) = pole(k)
        end do
      end do
    end if
    ! --------------------------------- FFSL -----------------------------------
    call adv_fill_vhalo(pt)
    call adv_calc_tracer_vflx(block%adv_batch_pt, pt, ptf_lev, dt)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dpt%d(i,j,k) = dpt%d(i,j,k) - (ptf_lev%d(i,j,k+1) - ptf_lev%d(i,j,k))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_grad_ptf')

  end subroutine calc_grad_ptf

  subroutine calc_dmgsdt(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    call perf_start('calc_dmgsdt')

    associate (mesh => block%mesh   , &
               dmf  => block%aux%dmf, & ! in
               dmgs => dtend%dmgs   )   ! out
    dmgs%d = 0
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dmgs%d(i,j) = dmgs%d(i,j) - dmf%d(i,j,k)
        end do
      end do
    end do
    end associate

    call perf_stop('calc_dmgsdt')

  end subroutine calc_dmgsdt

  subroutine calc_wedudlev_wedvdlev(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k

    ! Follow SB81 vertical advection discretization.

    call perf_start('calc_wedudlev_wedvdlev')

    associate (mesh       => block%mesh          , &
               u          => dstate%u_lon        , & ! in
               v          => dstate%v_lat        , & ! in
               dmg_lon    => block%aux%dmg_lon   , & ! in
               dmg_lat    => block%aux%dmg_lat   , & ! in
               we_lev_lon => block%aux%we_lev_lon, & ! in
               we_lev_lat => block%aux%we_lev_lat, & ! in
               du         => dtend%du            , & ! out
               dv         => dtend%dv            )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp = -(                                                &
            we_lev_lon%d(i,j,k+1) * (u%d(i,j,k+1) - u%d(i,j,k)) - &
            we_lev_lon%d(i,j,k  ) * (u%d(i,j,k-1) - u%d(i,j,k))   &
          ) / dmg_lon%d(i,j,k) / 2.0_r8
          du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_wedudeta%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp = -(                                                &
            we_lev_lat%d(i,j,k+1) * (v%d(i,j,k+1) - v%d(i,j,k)) - &
            we_lev_lat%d(i,j,k  ) * (v%d(i,j,k-1) - v%d(i,j,k))   &
          ) / dmg_lat%d(i,j,k) / 2.0_r8
          dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_wedvdeta%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_wedudlev_wedvdlev')

  end subroutine calc_wedudlev_wedvdlev

end module operators_mod
