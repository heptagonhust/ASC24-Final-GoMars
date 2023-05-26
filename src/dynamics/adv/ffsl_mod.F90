module ffsl_mod

  use const_mod
  use namelist_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use adv_batch_mod
  use ppm_mod
  use limiter_mod

  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_hflx
  public ffsl_calc_mass_vflx
  public ffsl_calc_tracer_hflx
  public ffsl_calc_tracer_vflx
  public ffsl_calc_tracer_vflx_lev

  interface
    subroutine hflx_interface(block, batch, u, v, mx, my, mfx, mfy)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: u  (block%mesh%half_ims:block%mesh%half_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
      real(r8), intent(in ) :: v  (block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%half_jms:block%mesh%half_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
      real(r8), intent(in ) :: mx (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%full_kms:block%filter_mesh%full_kme)
      real(r8), intent(in ) :: my (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%full_kms:block%filter_mesh%full_kme)
      real(r8), intent(out) :: mfx(block%mesh%half_ims:block%mesh%half_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
      real(r8), intent(out) :: mfy(block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%half_jms:block%mesh%half_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
    end subroutine hflx_interface
    subroutine vflx_interface(block, batch, w, m, mfz)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: w  (block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%half_kms:block%mesh%half_kme)
      real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%full_kms:block%filter_mesh%full_kme)
      real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%half_kms:block%mesh%half_kme)
    end subroutine vflx_interface
    subroutine vflx_lev_interface(block, batch, w, m, mfz)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: w  (block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
      real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%half_kms:block%filter_mesh%half_kme)
      real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
    end subroutine vflx_lev_interface
  end interface

  procedure(hflx_interface    ), pointer :: hflx     => null()
  procedure(vflx_interface    ), pointer :: vflx     => null()
  procedure(vflx_lev_interface), pointer :: vflx_lev => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      hflx     => hflx_van_leer
      vflx     => vflx_van_leer
    case ('ppm')
      hflx     => hflx_ppm
      vflx     => vflx_ppm
      vflx_lev => vflx_ppm_lev
    case default
      call log_error('Invalid ffsl_flux_type ' // trim(ffsl_flux_type) // '!', pid=proc%id)
    end select

    call limiter_init()

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_hflx(block, batch, m, mfx, mfy, dt)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfx(block%mesh%half_ims:block%mesh%half_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(out) :: mfy(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%half_jms:block%mesh%half_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)
    real(r8) dt_

    dt_ = merge(dt, batch%dt, present(dt))

    associate (mesh => block%mesh, &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx(block, batch, u, v, m, m, mfx, mfy)
    call fill_halo(block%halo, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(block%halo, mfy, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., east_halo=.false., north_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          ! Subtract divergence terms from flux to form advective operators.
          mx(i,j,k) = m(i,j,k) - 0.5_r8 * (          &
            (                                        &
              mfx(i,j,k) - mfx(i-1,j,k)              &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx(i,j,k) * m(i,j,k)                   &
          ) * dt_
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (     &
            (                                   &
              mfy(i,j  ,k) * mesh%le_lat(j  ) - &
              mfy(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j) -             &
            divy(i,j,k) * m(i,j,k)              &
          ) * dt_
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k)) * dt_
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k)) * dt_
        end do
      end do
    end if
    call fill_halo(block%filter_halo, mx, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false., cross_pole=.true.)
    call fill_halo(block%filter_halo, my, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false., cross_pole=.true.)
    ! Run outer flux form operators.
    call hflx(block, batch, u, v, my, mx, mfx, mfy)
    end associate

  end subroutine ffsl_calc_mass_hflx

  subroutine ffsl_calc_mass_vflx(block, batch, m, mfz, dt)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%half_kms:block%mesh%half_kme)
    real(r8), intent(in), optional :: dt

    call vflx(block, batch, batch%we, m, mfz)

  end subroutine ffsl_calc_mass_vflx

  subroutine ffsl_calc_tracer_hflx(block, batch, q, qmfx, qmfy, dt)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q    (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: qmfx (block%mesh%half_ims:block%mesh%half_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(out) :: qmfy (block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%half_jms:block%mesh%half_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)
    real(r8) dt_

    dt_ = merge(dt, batch%dt, present(dt))

    associate (mesh => block%mesh, &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx(block, batch, u, v, q, q, qmfx, qmfy)
    call fill_halo(block%halo, qmfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block%halo, qmfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          ! Subtract divergence terms from flux to form advective operators.
          qx(i,j,k) = q(i,j,k) - 0.5_r8 * (          &
            (                                        &
              qmfx(i,j,k) - qmfx(i-1,j,k)            &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx(i,j,k) * q(i,j,k)                   &
          ) * dt_
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (      &
            (                                    &
              qmfy(i,j  ,k) * mesh%le_lat(j  ) - &
              qmfy(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j) -              &
            divy(i,j,k) * q(i,j,k)               &
          ) * dt_
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = qmfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k)) * dt_
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = qmfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k)) * dt_
        end do
      end do
    end if
    call fill_halo(block%filter_halo, qx, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false., cross_pole=.true.)
    call fill_halo(block%filter_halo, qy, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false., cross_pole=.true.)
    ! Run outer flux form operators.
    call hflx(block, batch, mfx, mfy, qy, qx, qmfx, qmfy)
    end associate

  end subroutine ffsl_calc_tracer_hflx

  subroutine ffsl_calc_tracer_vflx(block, batch, q, qmfz, dt)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q   (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                  block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                  block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: qmfz(block%mesh%full_ims:block%mesh%full_ime, &
                                  block%mesh%full_jms:block%mesh%full_jme, &
                                  block%mesh%half_kms:block%mesh%half_kme)
    real(r8), intent(in), optional :: dt

    call vflx(block, batch, batch%we, q, qmfz)

  end subroutine ffsl_calc_tracer_vflx

  subroutine ffsl_calc_tracer_vflx_lev(block, batch, q, qmfz, dt)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q   (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                  block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                  block%filter_mesh%half_kms:block%filter_mesh%half_kme)
    real(r8), intent(out) :: qmfz(block%mesh%full_ims:block%mesh%full_ime, &
                                  block%mesh%full_jms:block%mesh%full_jme, &
                                  block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in), optional :: dt

    call vflx_lev(block, batch, batch%we, q, qmfz)

  end subroutine ffsl_calc_tracer_vflx_lev

  subroutine hflx_van_leer(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_ims:block%mesh%half_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in ) :: v  (block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%half_jms:block%mesh%half_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in ) :: mx (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(in ) :: my (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfx(block%mesh%half_ims:block%mesh%half_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(out) :: mfy(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%half_jms:block%mesh%half_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)

    integer i, j, k, iu, ju, ci
    real(r8) cf, dm

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    do k = mesh%full_kds, mesh%full_kde
      ! Along x-axis
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) + dm * 0.5_r8 * (1 - cf)) + sum(mx(i+1-ci:i,j,k))) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci + 1
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) - dm * 0.5_r8 * (1 + cf)) - sum(mx(i+1:i-ci,j,k))) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          cf = cfly(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          dm = slope(my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k))
          mfy(i,j,k) = v(i,j,k) * (my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    end associate

  end subroutine hflx_van_leer

  subroutine vflx_van_leer(block, batch, w, m, mfz)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: w  (block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%half_kms:block%mesh%half_kme)
    real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%half_kms:block%mesh%half_kme)

    integer i, j, k, ku, ci
    real(r8) cf, dm

    associate (mesh => block%mesh, &
               cflz => batch%cflz)   ! in
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ci = int(cflz(i,j,k))
          cf = cflz(i,j,k) - ci
          if (cflz(i,j,k) > 0) then
            ku = k - ci - 1
            dm = slope(m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1))
            mfz(i,j,k) = w(i,j,k) * (cf * (m(i,j,ku) + dm * 0.5_r8 * (1 - cf)) + sum(m(i,j,k-ci:k-1))) / cflz(i,j,k)
          else if (cflz(i,j,k) < 0) then
            ku = k - ci
            dm = slope(m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1))
            mfz(i,j,k) = w(i,j,k) * (cf * (m(i,j,ku) - dm * 0.5_r8 * (1 + cf)) - sum(m(i,j,k:k-ci-1))) / cflz(i,j,k)
          else
            mfz(i,j,k) = 0
          end if
        end do
      end do
    end do
    end associate

  end subroutine vflx_van_leer

  subroutine hflx_ppm(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_ims:block%mesh%half_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in ) :: v  (block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%half_jms:block%mesh%half_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in ) :: mx (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(in ) :: my (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfx(block%mesh%half_ims:block%mesh%half_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(out) :: mfy(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%half_jms:block%mesh%half_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)

    integer i, j, k, iu, ju, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    do k = mesh%full_kds, mesh%full_kde
      ! Along x-axis
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci
            call ppm(mx(iu-2,j,k), mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k), mx(iu+2,j,k), ml, dm, m6)
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) =  u(i,j,k) * (sum(mx(i+1-ci:i,j,k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci + 1
            call ppm(mx(iu-2,j,k), mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k), mx(iu+2,j,k), ml, dm, m6)
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) = -u(i,j,k) * (sum(mx(i+1:i-ci,j,k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          if (cfly(i,j,k) > 0) then
            ju = j
            call ppm(my(i,ju-2,k), my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k), my(i,ju+2,k), ml, dm, m6)
            s1 = 1 - cfly(i,j,k)
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) =  v(i,j,k) * (ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else if (cfly(i,j,k) < 0) then
            ju = j + 1
            call ppm(my(i,ju-2,k), my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k), my(i,ju+2,k), ml, dm, m6)
            s1 = 0
            s2 = -cfly(i,j,k)
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) = -v(i,j,k) * (ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else
            mfy(i,j,k) = 0
          end if
        end do
      end do
    end do
    end associate

  end subroutine hflx_ppm

  subroutine vflx_ppm(block, batch, w, m, mfz)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: w  (block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%half_kms:block%mesh%half_kme)
    real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)
    real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%half_kms:block%mesh%half_kme)

    integer i, j, k, ku, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6

    associate (mesh => block%mesh, &
               cflz => batch%cflz)   ! in
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ci = int(cflz(i,j,k))
          cf = cflz(i,j,k) - ci
          if (cflz(i,j,k) > 0) then
            ku = k - ci - 1
            call ppm(m(i,j,ku-2), m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1), m(i,j,ku+2), ml, dm, m6)
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfz(i,j,k) =  w(i,j,k) * (sum(m(i,j,k-ci:k-1)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz(i,j,k)
          else if (cflz(i,j,k) < 0) then
            ku = k - ci
            call ppm(m(i,j,ku-2), m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1), m(i,j,ku+2), ml, dm, m6)
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfz(i,j,k) = -w(i,j,k) * (sum(m(i,j,k:k-ci-1)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz(i,j,k)
          else
            mfz(i,j,k) = 0
          end if
        end do
      end do
    end do
    end associate

  end subroutine vflx_ppm

  subroutine vflx_ppm_lev(block, batch, w, m, mfz)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: w  (block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)
    real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%half_kms:block%filter_mesh%half_kme)
    real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                 block%mesh%full_jms:block%mesh%full_jme, &
                                 block%mesh%full_kms:block%mesh%full_kme)

    integer i, j, k, ku, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6

    associate (mesh => block%mesh, &
               cflz => batch%cflz)   ! in
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          ci = int(cflz(i,j,k))
          cf = cflz(i,j,k) - ci
          if (cflz(i,j,k) > 0) then
            ku = k - ci
            call ppm(m(i,j,ku-2), m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1), m(i,j,ku+2), ml, dm, m6)
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfz(i,j,k) =  w(i,j,k) * (sum(m(i,j,k-ci+1:k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz(i,j,k)
          else if (cflz(i,j,k) < 0) then
            ku = k - ci + 1
            call ppm(m(i,j,ku-2), m(i,j,ku-1), m(i,j,ku), m(i,j,ku+1), m(i,j,ku+2), ml, dm, m6)
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfz(i,j,k) = -w(i,j,k) * (sum(m(i,j,k+1:k-ci)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz(i,j,k)
          else
            mfz(i,j,k) = (m(i,j,k) + m(i,j,k+1)) * 0.5_r8
          end if
        end do
      end do
    end do
    end associate

  end subroutine vflx_ppm_lev

end module ffsl_mod
