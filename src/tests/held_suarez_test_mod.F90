module held_suarez_test_mod

  use flogger
  use const_mod, only: r8, Rd_o_cpd
  use namelist_mod
  use formula_mod
  use vert_coord_mod
  use block_mod
  use latlon_parallel_mod
  use rossby_haurwitz_wave_3d_test_mod

  implicit none

  private

  public held_suarez_test_set_ic
  public held_suarez_test_apply_forcing

  real(r8), parameter :: sig_b    = 0.7_r8
  real(r8), parameter :: kf       = 1.0_r8 / 86400.0_r8   ! s-1
  real(r8), parameter :: ka       = 0.025_r8 / 86400.0_r8 ! s-1
  real(r8), parameter :: ks       = 0.25_r8 / 86400.0_r8  ! s-1
  real(r8), parameter :: dt_lat   = 60.0_r8               ! K
  real(r8), parameter :: dpt_lev  = 10.0_r8               ! K
  real(r8), parameter :: p0       = 1.0e5_r8              ! Pa

contains

  subroutine held_suarez_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) random

    call rossby_haurwitz_wave_3d_test_set_ic(block)

    call random_seed()

    associate (mesh  => block%mesh         , &
               mgs   => block%dstate(1)%mgs, &
               mg    => block%dstate(1)%mg , &
               t     => block%dstate(1)%t  , &
               pt    => block%dstate(1)%pt )
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        call random_number(random)
        mgs%d(i,j) = mgs%d(i,j) - (0.5_r8 + random) * mesh%full_cos_lat(j)**2
      end do
    end do
    call fill_halo(mgs)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          call random_number(random)
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), 0.0_r8) - (0.5_r8 + random) * mesh%full_cos_lat(j)**2
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine held_suarez_test_set_ic

  subroutine held_suarez_test_apply_forcing(block, dt, dstate)

    type(block_type), intent(in) :: block
    real(r8), intent(in) :: dt
    type(dstate_type), intent(inout) :: dstate

    real(r8) kv, kt, teq, p_p0
    integer i, j, k

    associate (mesh  => block%mesh  , &
               u     => dstate%u_lon, &
               v     => dstate%v_lat, &
               mg    => dstate%mg   , &
               t     => dstate%t    , &
               pt    => dstate%pt   )
    do k = mesh%full_kds, mesh%full_kde
      kv = kf * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b))
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u%d(i,j,k) = u%d(i,j,k) - dt * kv * u%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(u)

    do k = mesh%full_kds, mesh%full_kde
      kv = kf * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b))
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v%d(i,j,k) = v%d(i,j,k) - dt * kv * v%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(v)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        kt = ka + (ks - ka) * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b)) * mesh%full_cos_lat(j)**4
        do i = mesh%full_ids, mesh%full_ide
          p_p0 = mg%d(i,j,k) / p0
          teq = max(200.0_r8, (315.0_r8 - dt_lat * mesh%full_sin_lat(j)**2 - dpt_lev * log(p_p0) * mesh%full_cos_lat(j)**2) * p_p0**Rd_o_cpd)
          pt%d(i,j,k) = pt%d(i,j,k) - dt * kt * pt%d(i,j,k) * (1.0_r8 - teq / t%d(i,j,k))
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine held_suarez_test_apply_forcing

end module held_suarez_test_mod
