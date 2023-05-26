module nh_mod

  use const_mod
  use namelist_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod
  use interp_mod
  use math_mod
  use debug_mod

  implicit none

  private

  public nh_prepare
  public nh_solve

contains

  subroutine nh_prepare(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call calc_rhod (blocks(iblk), blocks(iblk)%dstate(1))
      call calc_p    (blocks(iblk), blocks(iblk)%dstate(1))
      call interp_p  (blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine nh_prepare

  subroutine nh_solve(block, dtend, old_state, star_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(dtend_type), intent(inout) :: dtend
    type(dstate_type), intent(in) :: old_state
    type(dstate_type), intent(inout) :: star_state
    type(dstate_type), intent(inout) :: new_state
    real(r8), intent(in) :: dt

  end subroutine nh_solve

  subroutine calc_rhod(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    ! Diagnose dry air density from hydrostatic equation.
    associate (mesh   => block%mesh   , &
               gz_lev => dstate%gz_lev, & ! in
               dmg    => dstate%dmg   , & ! in
               rhod   => dstate%rhod  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          rhod(i,j,k) = - dmg(i,j,k) / (gz_lev(i,j,k+1) - gz_lev(i,j,k))
          if (rhod(i,j,k) <= 0) then
            print *, mesh%full_lon_deg(i), '(', to_str(i), ')', mesh%full_lat_deg(j), '(', to_str(j), ')', k
            print *, 'Model is instable!'
            call process_stop(1)
          end if
        end do
      end do
    end do
    end associate

  end subroutine calc_rhod

  subroutine calc_p(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    real(r8), parameter :: p0 = 1.0e5_r8
    integer i, j, k

    associate (mesh => block%mesh, & ! in
               rhod => dstate%rhod, & ! in
               pt   => dstate%pt  , & ! in
               p    => dstate%p)      ! out
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            p(i,j,k) = p0 * (Rd * pt(i,j,k) * rhod(i,j,k) / p0)**cpd_o_cvd
            if (debug_is_inf(p(i,j,k))) then
              print *, i, j, k, pt(i,j,k), rhod(i,j,k)
              stop 'NaN p!'
            end if
          end do
        end do
      end do
    end associate

  end subroutine calc_p

  subroutine calc_linearized_p(block, old_state, new_state)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(in) :: old_state
    type(dstate_type), intent(inout) :: new_state

    integer i, j, k

    associate (mesh       => block%mesh      , &
               old_dmg    => old_state%dmg   , & ! in
               new_dmg    => new_state%dmg   , & ! in
               old_p      => old_state%p     , & ! in
               new_p      => new_state%p     , & ! out
               old_gz_lev => old_state%gz_lev, & ! in
               new_gz_lev => new_state%gz_lev, & ! in
               old_pt     => old_state%pt    , & ! in
               new_pt     => new_state%pt)       ! in
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            new_p(i,j,k) = old_p(i,j,k) * (1.0_r8 + cpd_o_cvd * ( &
              (new_dmg(i,j,k) * new_pt(i,j,k)) /                  &
              (old_dmg(i,j,k) * old_pt(i,j,k)) -                  &
              (new_gz_lev(i,j,k+1) - new_gz_lev(i,j,k)) /         &
              (old_gz_lev(i,j,k+1) - old_gz_lev(i,j,k))           &
            ))
          end do
        end do
      end do
    end associate

  end subroutine calc_linearized_p

  subroutine interp_p(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    associate (mesh       => block%mesh      , &
               p          => dstate%p        , & ! in
               p_lev      => dstate%p_lev    , & ! out
               p_lev_lon  => dstate%p_lev_lon, & ! out
               p_lev_lat  => dstate%p_lev_lat)    ! out
      call interp_cell_to_lev_edge(mesh, p, p_lev)
      call fill_halo(block%halo, p_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      call interp_lev_edge_to_lev_lon_edge(mesh, p_lev, p_lev_lon)
      call interp_lev_edge_to_lev_lat_edge(mesh, p_lev, p_lev_lat)
      call fill_halo(block%halo, p_lev_lon, full_lon=.false., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false., north_halo=.false.)
    end associate

  end subroutine interp_p

  subroutine apply_bc_w_lev(block, star_state, new_state)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: star_state
    type(dstate_type), intent(inout) :: new_state

    real(r8) us_dzsdlon, vs_dzsdlat
    integer i, j, k

    associate (mesh           => block%mesh          , &
               u_lev_lon      => star_state%u_lev_lon, & ! in
               v_lev_lat      => star_state%v_lev_lat, & ! in
               dzsdlon        => block%static%dzsdlon, & ! in
               dzsdlat        => block%static%dzsdlat, & ! in
               w_lev          => new_state%w_lev)        ! out
    k = mesh%half_kde
    do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do i = mesh%full_ids, mesh%full_ide
        us_dzsdlon = (u_lev_lon(i-1,j,k) * dzsdlon(i-1,j) + &
                      u_lev_lon(i  ,j,k) * dzsdlon(i  ,j)) * 0.5_r8
        vs_dzsdlat = (v_lev_lat(i,j-1,k) * dzsdlat(i,j-1) + &
                      v_lev_lat(i,j  ,k) * dzsdlat(i,j  )) * 0.5_r8
        w_lev(i,j,k) = us_dzsdlon + vs_dzsdlat
      end do
    end do
    end associate

  end subroutine apply_bc_w_lev

  subroutine implicit_w_solver(block, dtend, old_state, star_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(dtend_type), intent(in) :: dtend
    type(dstate_type), intent(in) :: old_state
    type(dstate_type), intent(in) :: star_state
    type(dstate_type), intent(inout) :: new_state
    real(r8), intent(in) :: dt

    real(r8) w1 (block%mesh%half_kms:block%mesh%half_kme)
    real(r8) gz1(block%mesh%half_kms:block%mesh%half_kme)
    real(r8) dgz(block%mesh%full_kms:block%mesh%full_kme)
    real(r8) dp1, gdtbeta, gdt1mbeta, gdtbeta2gam
    real(r8) a(global_mesh%half_nlev)
    real(r8) b(global_mesh%half_nlev)
    real(r8) c(global_mesh%half_nlev)
    real(r8) d(global_mesh%half_nlev)
    integer i, j, k

    call apply_bc_w_lev(block, star_state, new_state)

    !
    ! ϕ¹ = ϕⁿ - Δt adv_ϕ* + g Δt (1 - β) w*
    !
    !
    ! w¹ = wⁿ - Δt adv_w* - g Δt + g Δt (1 - β) (∂p/∂π)*
    !
    ! Linearized dstate of ideal gas
    !
    ! 𝜹pⁿ⁺¹ ≈ 𝜹pⁿ + 𝜹(𝜸 pⁿ (𝜹𝜋 θ)ⁿ⁺¹ / (𝜹𝜋 θ)ⁿ) - 𝜹(𝜸 pⁿ 𝜹ϕ¹ / 𝜹ϕⁿ) - 𝜹(𝜸 pⁿ g Δt β 𝜹wⁿ⁺¹ / 𝜹φⁿ)
    !         -----------------------------------------------------
    !                                dp1
    !
    associate (mesh        => block%mesh        , &
               beta        => implicit_w_wgt    , &
               adv_gz      => dtend%adv_gz      , & ! FIXME: After test success, merge advection tends togethor.
               adv_w       => dtend%adv_w       , & !
               old_p       => old_state%p       , &
               star_p      => star_state%p      , &
               star_p_lev  => star_state%p_lev  , &
               old_w_lev   => old_state%w_lev   , &
               star_w_lev  => star_state%w_lev  , &
               new_w_lev   => new_state%w_lev   , &
               star_m_lev  => star_state%dmg_lev, &
               new_dmg_lev => new_state%dmg_lev , &
               old_gz_lev  => old_state%gz_lev  , &
               star_gz_lev => star_state%gz_lev , &
               new_gz_lev  => new_state%gz_lev  , &
               old_dmg     => old_state%dmg     , &
               new_dmg     => new_state%dmg     , &
               old_pt      => old_state%pt      , &
               new_pt      => new_state%pt)
      ! last: n, old: *, new: n + 1
      gdtbeta     = g * dt * beta
      gdt1mbeta   = g * dt * (1 - beta)
      gdtbeta2gam = (g * dt * beta)**2 * cpd_o_cvd
      ! FIXME: Two Poles may skip the duplicate calculation?
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          do k = mesh%full_kds, mesh%full_kde
            dgz(k) = old_gz_lev(i,j,k+1) - old_gz_lev(i,j,k)
          end do
          gz1 = old_gz_lev(i,j,:) - dt * adv_gz(i,j,:) + gdt1mbeta * star_w_lev(i,j,:)
          w1  = old_w_lev (i,j,:) - dt * adv_w (i,j,:) - g * dt
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            w1(k) = w1(k) + gdt1mbeta * (star_p(i,j,k) - star_p(i,j,k-1)) / star_m_lev(i,j,k)
          end do
          ! Top boundary
          k = mesh%half_kds
          w1(k) = w1(k) + gdt1mbeta * (star_p(i,j,k) - star_p_lev(i,j,k)) / star_m_lev(i,j,k)
          ! Bottom boundary
          k = mesh%half_kde
          w1(k) = w1(k) + gdt1mbeta * (star_p_lev(i,j,k) - star_p(i,j,k-1)) / star_m_lev(i,j,k)
          ! Use linearized dstate of ideal gas to calculate the first part of ∂pⁿ⁺¹ (i.e. dp1).
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            dp1 = (old_p(i,j,k) - old_p(i,j,k-1)) + cpd_o_cvd * ((                                       &
              old_p(i,j,k  ) * new_dmg(i,j,k  ) * new_pt(i,j,k  ) / old_dmg(i,j,k  ) / old_pt(i,j,k  ) - &
              old_p(i,j,k-1) * new_dmg(i,j,k-1) * new_pt(i,j,k-1) / old_dmg(i,j,k-1) / old_pt(i,j,k-1)   &
            ) - (                                                                                        &
              old_p(i,j,k  ) * (gz1(k+1) - gz1(k  )) / dgz(k  ) -                                        &
              old_p(i,j,k-1) * (gz1(k  ) - gz1(k-1)) / dgz(k-1)                                          &
            ))
            w1(k) = w1(k) + gdtbeta * dp1 / new_dmg_lev(i,j,k)
          end do
          ! Set coefficients for implicit solver.
          a(1) = 0.0_r8
          b(1) = 1.0_r8
          c(1) = 0.0_r8
          d(1) = 0.0_r8 ! Top w is set to zero.
          do k = mesh%half_kds + 1, mesh%half_kde - 1
            a(k) = gdtbeta2gam * old_p(i,j,k-1) / dgz(k-1)
            b(k) = new_dmg_lev(i,j,k) - gdtbeta2gam * (old_p(i,j,k) / dgz(k) + old_p(i,j,k-1) / dgz(k-1))
            c(k) = gdtbeta2gam * old_p(i,j,k  ) / dgz(k  )
            d(k) = new_dmg_lev(i,j,k) * w1(k)
          end do
          a(mesh%half_nlev) = 0.0_r8
          b(mesh%half_nlev) = 1.0_r8
          c(mesh%half_nlev) = 0.0_r8
          d(mesh%half_nlev) = new_w_lev(i,j,mesh%half_nlev)
          call tridiag_thomas(a, b, c, d, new_w_lev(i,j,:))

          call rayleigh_damp_w(dt, star_gz_lev(i,j,:), new_w_lev(i,j,:))

          ! Update gz after w is solved.
          do k = mesh%half_kds, mesh%half_kde - 1
            new_gz_lev(i,j,k) = gz1(k) + gdtbeta * new_w_lev(i,j,k)
          end do
        end do
      end do
      call fill_halo(block%halo, new_w_lev , full_lon=.true., full_lat=.true., full_lev=.false.)
      call fill_halo(block%halo, new_gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
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

end module nh_mod
