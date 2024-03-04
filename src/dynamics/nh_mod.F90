! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module nh_mod

  use const_mod
  use namelist_mod
  use block_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use latlon_parallel_mod
  use process_mod
  use perf_mod
  use interp_mod
  use math_mod
  use adv_mod
  use tracer_mod
  use operators_mod
  use filter_mod

  implicit none

  private

  public nh_solve

contains

  subroutine nh_solve(block, old_dstate, star_dstate, new_dstate, dt)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(inout) :: star_dstate
    type(dstate_type), intent(inout) :: new_dstate
    real(r8), intent(in) :: dt

    call interp_wind(block, star_dstate)
    call calc_adv_lev(block, star_dstate%w_lev , block%aux%adv_w_lev , star_dstate%dmg_lev, dt)
    call calc_adv_lev(block, star_dstate%gz_lev, block%aux%adv_gz_lev, star_dstate%dmg_lev, dt)
    call implicit_w_solver(block, old_dstate, star_dstate, new_dstate, dt)
    call calc_p(block, old_dstate, new_dstate, dt)
    call average_run(new_dstate%gz_lev, new_dstate%gz)
    call fill_halo(new_dstate%gz, west_halo=.false., south_halo=.false.)

  end subroutine nh_solve

  subroutine interp_wind(block, dstate)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    associate (u_lon       => dstate%u_lon         , & ! in
               v_lat       => dstate%v_lat         , & ! in
               we_lev      => dstate%we_lev        , & ! in
               dmg         => dstate%dmg           , & ! in
               dmg_lev     => dstate%dmg_lev       , & ! in
               u_lev_lon   => block%aux%u_lev_lon  , & ! out
               v_lev_lat   => block%aux%v_lev_lat  , & ! out
               we          => dstate%we            , & ! out
               mfx_lev_lon => block%aux%mfx_lev_lon, & ! out
               mfy_lev_lat => block%aux%mfy_lev_lat, & ! out
               dmf_lev     => block%aux%dmf_lev    )   ! out
    call interp_run(u_lon, u_lev_lon)
    call interp_run(v_lat, v_lev_lat)
    call interp_run(dmg_lev, mfx_lev_lon)
    mfx_lev_lon%d = mfx_lev_lon%d * u_lev_lon%d
    call interp_run(dmg_lev, mfy_lev_lat)
    mfy_lev_lat%d = mfy_lev_lat%d * v_lev_lat%d
    call interp_run(we_lev, we)
    call fill_halo(u_lev_lon)
    call fill_halo(v_lev_lat)
    call fill_halo(mfx_lev_lon)
    call fill_halo(mfy_lev_lat)
    call block%adv_batch_nh%set_wind(u_lev_lon, v_lev_lat, we, mfx_lev_lon, mfy_lev_lat, dmg)
    call div_operator(mfx_lev_lon, mfy_lev_lat, dmf_lev)
    end associate

  end subroutine interp_wind

  subroutine calc_adv_lev(block, q_lev, dqdt_lev, dmg_lev, dt)

    type(block_type         ), intent(inout) :: block
    type(latlon_field3d_type), intent(inout) :: q_lev
    type(latlon_field3d_type), intent(inout) :: dqdt_lev
    type(latlon_field3d_type), intent(in   ) :: dmg_lev
    real(r8), intent(in) :: dt

    integer i, j, k

    associate (mesh        => block%mesh                , &
               dmf_lev     => block%aux%dmf_lev         , & ! in
               we          => block%adv_batch_nh%we     , & ! in
               qmf_lon     => block%adv_batch_nh%qmf_lon, &
               qmf_lat     => block%adv_batch_nh%qmf_lat, &
               qmf_lev     => block%adv_batch_nh%qmf_lev)
    call adv_calc_tracer_hflx(block%adv_batch_nh, q_lev, qmf_lon, qmf_lat, dt)
    call fill_halo(qmf_lon, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(qmf_lat, west_halo=.false., east_halo=.false., north_halo=.false.)
    call div_operator(qmf_lon, qmf_lat, dqdt_lev)
    ! Remove horizontal mass flux divergence part.
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dqdt_lev%d(i,j,k) = -dqdt_lev%d(i,j,k) + q_lev%d(i,j,k) * dmf_lev%d(i,j,k)
        end do
      end do
    end do
    call adv_fill_vhalo(q_lev)
    call adv_calc_tracer_vflx(block%adv_batch_nh, q_lev, qmf_lev, dt)
    ! Remove vertical mass flux divergence part.
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dqdt_lev%d(i,j,k) = dqdt_lev%d(i,j,k) - (       &
            qmf_lev%d(i,j,k) - qmf_lev%d(i,j,k-1) -       &
            q_lev%d(i,j,k) * (we%d(i,j,k) - we%d(i,j,k-1)))
        end do
      end do
    end do
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dqdt_lev%d(i,j,k) = dqdt_lev%d(i,j,k) / dmg_lev%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(dqdt_lev, south_halo=.false., north_halo=.false.)
    call filter_run(block%big_filter, dqdt_lev)
    end associate

  end subroutine calc_adv_lev

  subroutine apply_bc_w_lev(block, dstate)

    type(block_type ), intent(in   ) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    associate (mesh       => block%mesh          , &
               adv_gz_lev => block%aux%adv_gz_lev, & ! in
               w_lev      => dstate%w_lev        )   ! out
    k = mesh%half_kde
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        w_lev%d(i,j,k) = -adv_gz_lev%d(i,j,k) / g
      end do
    end do
    end associate

  end subroutine apply_bc_w_lev

  subroutine implicit_w_solver(block, old_dstate, star_dstate, new_dstate, dt)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(in   ) :: star_dstate
    type(dstate_type), intent(inout) :: new_dstate
    real(r8), intent(in) :: dt

    real(r8) w1 (block%mesh%half_kms:block%mesh%half_kme)
    real(r8) gz1(block%mesh%half_kms:block%mesh%half_kme)
    real(r8) dgz(block%mesh%full_kms:block%mesh%full_kme)
    real(r8) L  (block%mesh%half_kms:block%mesh%half_kme)
    real(r8) dp1, gdtbeta, gdt1mbeta, gdtbeta2gam
    real(r8) a(global_mesh%half_nlev)
    real(r8) b(global_mesh%half_nlev)
    real(r8) c(global_mesh%half_nlev)
    real(r8) d(global_mesh%half_nlev)
    integer i, j, k

    call apply_bc_w_lev(block, new_dstate)

    associate (mesh         => block%mesh              , &
               beta         => implicit_w_wgt          , &
               adv_gz_lev   => block%aux%adv_gz_lev    , & ! in
               adv_w_lev    => block%aux%adv_w_lev     , & ! in
               old_p        => old_dstate%p            , & ! in
               star_p       => star_dstate%p           , & ! in
               star_p_lev   => star_dstate%p_lev       , & ! in
               old_w_lev    => old_dstate%w_lev        , & ! in
               star_w_lev   => star_dstate%w_lev       , & ! in
               new_w_lev    => new_dstate%w_lev        , & ! out
               star_dmg_lev => star_dstate%dmg_lev     , & ! in
               new_dmg_lev  => new_dstate%dmg_lev      , & ! in
               old_gz_lev   => old_dstate%gz_lev       , & ! in
               star_gz_lev  => star_dstate%gz_lev      , & ! in
               new_gz_lev   => new_dstate%gz_lev       , & ! out
               old_dmg      => old_dstate%dmg          , & ! in
               new_dmg      => new_dstate%dmg          , & ! in
               old_pt       => old_dstate%pt           , & ! in
               new_pt       => new_dstate%pt           , & ! in
               qm_lev       => tracers(block%id)%qm_lev)   ! in
    ! last: n, old: *, new: n + 1
    gdtbeta     = g * dt * beta
    gdt1mbeta   = g * dt * (1 - beta)
    gdtbeta2gam = (g * dt * beta)**2 * cpd_o_cvd
    ! FIXME: Two Poles may skip the duplicate calculation?
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        do k = mesh%half_kds, mesh%half_kde
          L(k) = 1.0_r8 / (1 + qm_lev%d(i,j,k))
        end do
        do k = mesh%full_kds, mesh%full_kde
          dgz(k) = old_gz_lev%d(i,j,k+1) - old_gz_lev%d(i,j,k)
        end do
        !
        ! ϕ¹ = ϕⁿ + Δt adv_ϕ* + g Δt (1 - β) w*
        !
        gz1 = old_gz_lev%d(i,j,:) + dt * adv_gz_lev%d(i,j,:) + gdt1mbeta * star_w_lev%d(i,j,:)
        !
        ! w¹ = wⁿ + Δt adv_w* - g Δt + L g Δt (1 - β) 𝜹p* / 𝜹𝜋*
        !
        w1  = old_w_lev %d(i,j,:) + dt * adv_w_lev %d(i,j,:) - g * dt
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          w1(k) = w1(k) + gdt1mbeta * L(k) * (star_p%d(i,j,k) - star_p%d(i,j,k-1)) / star_dmg_lev%d(i,j,k)
        end do
        ! Top boundary
        k = mesh%half_kds
        w1(k) = w1(k) + gdt1mbeta * L(k) * (star_p%d(i,j,k) - star_p_lev%d(i,j,k)) / star_dmg_lev%d(i,j,k)
        ! Bottom boundary
        k = mesh%half_kde
        w1(k) = w1(k) + gdt1mbeta * L(k) * (star_p_lev%d(i,j,k) - star_p%d(i,j,k-1)) / star_dmg_lev%d(i,j,k)
        ! Use linearized dstate of ideal gas to calculate the first part of 𝜹pⁿ⁺¹ (i.e. dp1).
        !
        ! 𝜹pⁿ⁺¹ ≈ 𝜹pⁿ + 𝜸 𝜹(pⁿ (𝜹𝜋 θₘ)ⁿ⁺¹ / (𝜹𝜋 θₘ)ⁿ) - 𝜸 𝜹(pⁿ 𝜹ϕ¹ / 𝜹ϕⁿ) - β Δt g 𝜸 𝜹(pⁿ 𝜹wⁿ⁺¹ / 𝜹φⁿ)
        !         -----------------------------------------------------
        !                                dp1
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          dp1 = (old_p%d(i,j,k) - old_p%d(i,j,k-1)) + cpd_o_cvd * ((                                             &
            old_p%d(i,j,k  ) * new_dmg%d(i,j,k  ) * new_pt%d(i,j,k  ) / old_dmg%d(i,j,k  ) / old_pt%d(i,j,k  ) - &
            old_p%d(i,j,k-1) * new_dmg%d(i,j,k-1) * new_pt%d(i,j,k-1) / old_dmg%d(i,j,k-1) / old_pt%d(i,j,k-1)   &
          ) - (                                                                                                  &
            old_p%d(i,j,k  ) * (gz1(k+1) - gz1(k  )) / dgz(k  ) -                                                &
            old_p%d(i,j,k-1) * (gz1(k  ) - gz1(k-1)) / dgz(k-1)                                                  &
          ))
          w1(k) = w1(k) + gdtbeta * L(k) * dp1 / new_dmg_lev%d(i,j,k)
        end do

        ! Set coefficients for implicit solver.
        a(1) = 0.0_r8
        b(1) = 1.0_r8
        c(1) = 0.0_r8
        d(1) = 0.0_r8 ! Top w is set to zero.
        do k = mesh%half_kds + 1, mesh%half_kde - 1
          a(k) = gdtbeta2gam * L(k) * old_p%d(i,j,k-1) / dgz(k-1)
          c(k) = gdtbeta2gam * L(k) * old_p%d(i,j,k  ) / dgz(k  )
          b(k) = new_dmg_lev%d(i,j,k) - a(k) - c(k)
          d(k) = new_dmg_lev%d(i,j,k) * w1(k)
        end do
        a(mesh%half_nlev) = 0.0_r8
        b(mesh%half_nlev) = 1.0_r8
        c(mesh%half_nlev) = 0.0_r8
        d(mesh%half_nlev) = new_w_lev%d(i,j,mesh%half_nlev)
        call tridiag_thomas(a, b, c, d, new_w_lev%d(i,j,mesh%half_kds:mesh%half_kde))

        call rayleigh_damp_w(dt, star_gz_lev%d(i,j,mesh%half_kds:mesh%half_kde), &
                                   new_w_lev%d(i,j,mesh%half_kds:mesh%half_kde))

        ! Update gz after w is solved.
        do k = mesh%half_kds, mesh%half_kde - 1
          new_gz_lev%d(i,j,k) = gz1(k) + gdtbeta * new_w_lev%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(new_w_lev )
    call fill_halo(new_gz_lev)
    end associate

  end subroutine implicit_w_solver

  subroutine rayleigh_damp_w(dt, gz, w)

    real(r8), intent(in   ) :: dt
    real(r8), intent(in   ) :: gz(:)
    real(r8), intent(inout) :: w (:)

    real(r8) gzd, c
    integer k

    gzd = rayleigh_damp_top * g
    do k = 2, size(w) - 1
      if (gz(k) > gz(1) - gzd) then
        c = rayleigh_damp_w_coef * sin(pi05 * (1 - (gz(1) - gz(k)) / gzd))**2
        w(k) = w(k) / (1 + c * dt)
      else
        return
      end if
    end do

  end subroutine rayleigh_damp_w

  subroutine calc_p(block, old_dstate, new_dstate, dt)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(in) :: old_dstate
    type(dstate_type), intent(inout) :: new_dstate
    real(r8), intent(in) :: dt

    integer i, j, k

    call perf_start('calc_p')

    associate (mesh       => block%mesh              , & ! in
               old_p      => old_dstate%p            , & ! in
               new_p      => new_dstate%p            , & ! out
               new_p_lev  => new_dstate%p_lev        , & ! out
               old_pt     => old_dstate%pt           , & ! in
               new_pt     => new_dstate%pt           , & ! in
               old_dmg    => old_dstate%dmg          , & ! in
               new_dmg    => new_dstate%dmg          , & ! in
               old_gz_lev => old_dstate%gz_lev       , & ! in
               new_gz_lev => new_dstate%gz_lev       , & ! in
               old_w_lev  => old_dstate%w_lev        , & ! in
               new_w_lev  => new_dstate%w_lev        , & ! in
               adv_w_lev  => block%aux%adv_w_lev     , & ! in
               qm_lev     => tracers(block%id)%qm_lev)   ! in
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          new_p%d(i,j,k) = old_p%d(i,j,k) + cpd_o_cvd * old_p%d(i,j,k) * (                                &
            (new_dmg%d(i,j,k) * new_pt%d(i,j,k)) / (old_dmg%d(i,j,k) * old_pt%d(i,j,k)) -                 &
            (new_gz_lev%d(i,j,k) - new_gz_lev%d(i,j,k+1)) / (old_gz_lev%d(i,j,k) - old_gz_lev%d(i,j,k+1)) &
          )
          ! Do 3D divergence damping?
          new_p%d(i,j,k) = new_p%d(i,j,k) + 0.12_r8 * (new_p%d(i,j,k) - old_p%d(i,j,k))
        end do
      end do
    end do
    ! Calculate half level pressure from w equation.
    new_p_lev%d(:,:,1) = ptop
    do k = mesh%half_kds + 1, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          new_p_lev%d(i,j,k) = new_p_lev%d(i,j,k-1) + (            &
            ((new_w_lev%d(i,j,k-1) - old_w_lev%d(i,j,k-1)) / dt -  &
             adv_w_lev%d(i,j,k-1) + g) * (1 + qm_lev%d(i,j,k-1)) + &
            ((new_w_lev%d(i,j,k  ) - old_w_lev%d(i,j,k  )) / dt -  &
             adv_w_lev%d(i,j,k  ) + g) * (1 + qm_lev%d(i,j,k  ))   &
          ) * 0.5_r8 / g * new_dmg%d(i,j,k-1)
        end do
      end do
    end do
    call fill_halo(new_p_lev)
    end associate

    call perf_stop('calc_p')

  end subroutine calc_p

end module nh_mod
