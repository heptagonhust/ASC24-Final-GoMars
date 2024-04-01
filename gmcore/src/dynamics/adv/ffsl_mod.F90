! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module implements FFSL (Flux-Form Semi-Lagrangian) advection scheme on
!   the lat-lon grid.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module ffsl_mod

  use const_mod
  use namelist_mod
  use latlon_mesh_mod, only: global_mesh
  use latlon_field_types_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use adv_batch_mod
  use ppm_mod
  use limiter_mod
  use perf_mod
  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_hflx
  public ffsl_calc_mass_vflx
  public ffsl_calc_tracer_hflx
  public ffsl_calc_tracer_vflx

  interface
    subroutine hflx_interface(batch, u, v, mx, my, mfx, mfy)
      import  adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: u
      type(latlon_field3d_type), intent(in   ) :: v
      type(latlon_field3d_type), intent(in   ) :: mx
      type(latlon_field3d_type), intent(in   ) :: my
      type(latlon_field3d_type), intent(inout) :: mfx
      type(latlon_field3d_type), intent(inout) :: mfy
    end subroutine hflx_interface
    subroutine vflx_interface(batch, w, m, mfz)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: w
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfz
    end subroutine vflx_interface
  end interface

  procedure(hflx_interface), pointer :: hflx => null()
  procedure(vflx_interface), pointer :: vflx => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      hflx     => hflx_van_leer
      vflx     => vflx_van_leer
    case ('ppm')
      hflx     => hflx_ppm
      vflx     => vflx_ppm
    case default
      call log_error('Invalid ffsl_flux_type ' // trim(ffsl_flux_type) // '!', pid=proc%id)
    end select

    call limiter_init()

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_hflx(batch, m, mfx, mfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy
    real(r8), intent(in), optional :: dt

    integer i, j, k
    real(r8) work(m%mesh%full_ids:m%mesh%full_ide,m%mesh%full_nlev)
    real(r8) pole(m%mesh%full_nlev)
    real(r8) dt_opt

    call perf_start('ffsl_calc_mass_hflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh => m%mesh    , &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy  )   ! work array
    ! Run inner advective operators.
    call hflx(batch, u, v, m, m, mfx, mfy)
    call fill_halo(mfx, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(mfy, west_halo=.false., east_halo=.false., north_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          ! Subtract divergence terms from flux to form advective operators.
          mx%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (      &
            (                                        &
              mfx%d(i,j,k) - mfx%d(i-1,j,k)          &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx%d(i,j,k) * m%d(i,j,k)               &
          ) * dt_opt
          my%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (     &
            (                                       &
              mfy%d(i,j  ,k) * mesh%le_lat(j  ) -   &
              mfy%d(i,j-1,k) * mesh%le_lat(j-1)     &
            ) / mesh%area_cell(j) -                 &
            divy%d(i,j,k) * m%d(i,j,k)              &
          ) * dt_opt
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy%d(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = m%d(i,j,k)
          my%d(i,j,k) = m%d(i,j,k) - 0.5_r8 * (pole(k) - divy%d(i,j,k) * m%d(i,j,k)) * dt_opt
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = mfy%d(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          mx%d(i,j,k) = m%d(i,j,k)
          my%d(i,j,k) = m%d(i,j,k) + 0.5_r8 * (pole(k) - divy%d(i,j,k) * m%d(i,j,k)) * dt_opt
        end do
      end do
    end if
    call fill_halo(mx, west_halo=.false., east_halo=.false.)
    call fill_halo(my, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx(batch, u, v, my, mx, mfx, mfy)
    end associate

    call perf_stop('ffsl_calc_mass_hflx')

  end subroutine ffsl_calc_mass_hflx

  subroutine ffsl_calc_mass_vflx(batch, m, mfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz
    real(r8), intent(in), optional :: dt
    call perf_start('ffsl_calc_mass_vflx')
    call vflx(batch, batch%we, m, mfz)
    call perf_stop('ffsl_calc_mass_vflx')
  end subroutine ffsl_calc_mass_vflx

  subroutine ffsl_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfx
    type(latlon_field3d_type), intent(inout) :: qmfy
    real(r8), intent(in), optional :: dt

    integer ks, ke, i, j, k
    real(r8) work(q%mesh%full_ids:q%mesh%full_ide,q%mesh%half_nlev)
    real(r8) pole(q%mesh%half_nlev)
    real(r8) dt_opt

    call perf_start('ffsl_calc_tracer_hflx')

    dt_opt = batch%dt; if (present(dt)) dt_opt = dt

    associate (mesh => q%mesh    , &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy  , & ! work array
               q    => q         , & ! in
               qmfx => qmfx      , & ! out
               qmfy => qmfy      )   ! out
    ! Run inner advective operators.
    call hflx(batch, u, v, q, q, qmfx, qmfy)
    call fill_halo(qmfx, south_halo=.false., north_halo=.false., east_halo=.false.)
    call fill_halo(qmfy, west_halo=.false., east_halo=.false., north_halo=.false.)
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      ! Calculate intermediate tracer density due to advective operators.
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide
            ! Subtract divergence terms from flux to form advective operators.
            qx%d(i,j,k) = q%d(i,j,k) - 0.5_r8 * (      &
              (                                        &
                qmfx%d(i,j,k) - qmfx%d(i-1,j,k)        &
              ) * mesh%le_lon(j) / mesh%area_cell(j) - &
              divx%d(i,j,k) * q%d(i,j,k)               &
            ) * dt_opt
            qy%d(i,j,k) = q%d(i,j,k) - 0.5_r8 * (    &
              (                                      &
                qmfy%d(i,j  ,k) * mesh%le_lat(j  ) - &
                qmfy%d(i,j-1,k) * mesh%le_lat(j-1)   &
              ) / mesh%area_cell(j) -                &
              divy%d(i,j,k) * q%d(i,j,k)             &
            ) * dt_opt
          end do
        end do
      end do
      ! Handle the Pole boundary conditions.
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = qmfy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = q%d(i,j,k)
            qy%d(i,j,k) = q%d(i,j,k) - 0.5_r8 * (pole(k) - divy%d(i,j,k) * q%d(i,j,k)) * dt_opt
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = qmfy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j)
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            qx%d(i,j,k) = q%d(i,j,k)
            qy%d(i,j,k) = q%d(i,j,k) + 0.5_r8 * (pole(k) - divy%d(i,j,k) * q%d(i,j,k)) * dt_opt
          end do
        end do
      end if
    case ('vtx')
    end select
    call fill_halo(qx, west_halo=.false., east_halo=.false.)
    call fill_halo(qy, south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx(batch, mfx, mfy, qy, qx, qmfx, qmfy)
    end associate
    call perf_stop('ffsl_calc_tracer_hflx')
  end subroutine ffsl_calc_tracer_hflx

  subroutine ffsl_calc_tracer_vflx(batch, q, qmfz, dt)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: q
    type(latlon_field3d_type), intent(inout) :: qmfz
    real(r8), intent(in), optional :: dt
    call perf_start('ffsl_calc_tracer_vflx')
    call vflx(batch, batch%we, q, qmfz)
    call perf_stop('ffsl_calc_tracer_vflx')
  end subroutine ffsl_calc_tracer_vflx

  subroutine hflx_van_leer(batch, u, v, mx, my, mfx, mfy)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(in   ) :: mx
    type(latlon_field3d_type), intent(in   ) :: my
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy

    integer ks, ke, i, j, k, iu, ju, ci
    real(r8) cf, dm

    call perf_start('hflx_van_leer')
    associate (mesh => u%mesh    , &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')

      do k = ks, ke
        ! Along x-axis
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            ci = int(cflx%d(i,j,k))
            cf = cflx%d(i,j,k) - ci
            if (abs(cflx%d(i,j,k)) < 1.0e-16_r8) then
              mfx%d(i,j,k) = 0
            else if (cflx%d(i,j,k) > 0) then
              iu = i - ci
              dm = slope(mx%d(iu-1,j,k), mx%d(iu,j,k), mx%d(iu+1,j,k))
              mfx%d(i,j,k) = u%d(i,j,k) * (cf * (mx%d(iu,j,k) + dm * 0.5_r8 * (1 - cf)) + sum(mx%d(iu+1:i,j,k))) / cflx%d(i,j,k)
            else
              iu = i - ci + 1
              dm = slope(mx%d(iu-1,j,k), mx%d(iu,j,k), mx%d(iu+1,j,k))
              mfx%d(i,j,k) = u%d(i,j,k) * (cf * (mx%d(iu,j,k) - dm * 0.5_r8 * (1 + cf)) - sum(mx%d(i+1:iu-1,j,k))) / cflx%d(i,j,k)
            end if
          end do
        end do
        ! Along y-axis
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            cf = cfly%d(i,j,k)
            ju = merge(j, j + 1, cf > 0)
            dm = slope(my%d(i,ju-1,k), my%d(i,ju,k), my%d(i,ju+1,k))
            mfy%d(i,j,k) = v%d(i,j,k) * (my%d(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
          end do
        end do
      end do
    case ('vtx')
    end select
    end associate

    call perf_stop('hflx_van_leer')

  end subroutine hflx_van_leer

  subroutine vflx_van_leer(batch, w, m, mfz)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: w
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz

    integer i, j, k, ku, ci
    real(r8) cf, dm


    call perf_start('vflx_van_leer')

    associate (mesh => m%mesh    , &
               cflz => batch%cflz)   ! in
    select case (batch%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci - 1
              dm = slope(m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1))
              mfz%d(i,j,k) = w%d(i,j,k) * (cf * (m%d(i,j,ku) + dm * 0.5_r8 * (1 - cf)) + sum(m%d(i,j,ku+1:k-1))) / cflz%d(i,j,k)
            else
              ku = k - ci
              dm = slope(m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1))
              mfz%d(i,j,k) = w%d(i,j,k) * (cf * (m%d(i,j,ku) - dm * 0.5_r8 * (1 + cf)) - sum(m%d(i,j,k:ku-1))) / cflz%d(i,j,k)
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci
              dm = slope(m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1))
              mfz%d(i,j,k) = w%d(i,j,k) * (cf * (m%d(i,j,ku) + dm * 0.5_r8 * (1 - cf)) + sum(m%d(i,j,ku+1:k))) / cflz%d(i,j,k)
            else
              ku = k - ci + 1
              dm = slope(m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1))
              mfz%d(i,j,k) = w%d(i,j,k) * (cf * (m%d(i,j,ku) - dm * 0.5_r8 * (1 + cf)) - sum(m%d(i,j,k+1:ku-1))) / cflz%d(i,j,k)
            end if
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('vflx_van_leer')

  end subroutine vflx_van_leer
subroutine hflx_ppm(batch, u, v, mx, my, mfx, mfy)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(in   ) :: mx
    type(latlon_field3d_type), intent(in   ) :: my
    type(latlon_field3d_type), intent(inout) :: mfx
    type(latlon_field3d_type), intent(inout) :: mfy

    integer ks, ke, i, j, k, iu, ju, ci
    integer is, ie, js, je
    real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6
    real(8), pointer, dimension(:,:,:) :: cflx_d, mfx_d, mx_d, u_d, cfly_d, mfy_d, my_d, v_d

    ! interface
    !   subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)

    !   real(8), intent(in ) :: fm2
    !   real(8), intent(in ) :: fm1
    !   real(8), intent(in ) :: f
    !   real(8), intent(in ) :: fp1
    !   real(8), intent(in ) :: fp2
    !   real(8), intent(out) :: fl
    !   real(8), intent(out) :: df
    !   real(8), intent(out) :: f6
    !   !$omp declare target

    !   end subroutine ppm
    ! end interface
    call perf_start('hflx_ppm')
    associate (mesh => u%mesh    , &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    select case (batch%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole
      is = mesh%half_ids
      ie = mesh%half_ide
      cflx_d =>  cflx%d
      mfx_d  =>  mfx%d
      mx_d   =>  mx%d
      u_d    =>  u%d
      cfly_d =>  cfly%d
      mfy_d  =>  mfy%d
      my_d   =>  my%d
      v_d    =>  v%d
      !$omp target update to(cfly_d, my_d, v_d)
      !$omp target update to(cflx_d, mx_d, u_d)
      ! call cflx%HostToDevice()
      ! call mx%HostToDevice()
      ! call u%HostToDevice()
      ! call cfly%HostToDevice()
      ! call my%HostToDevice()
      ! call v%HostToDevice()
      !$omp target teams distribute parallel do private(i, j, k, iu, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(3)
      do k = ks, ke
        ! Along x-axis
        do j = js, je
          do i = is, ie
            ci = int(cflx_d(i,j,k))
            cf = cflx_d(i,j,k) - ci
            if (abs(cflx_d(i,j,k)) < 1.0e-16_r8) then
              mfx_d(i,j,k) = 0
            else if (cflx_d(i,j,k) > 0) then
              iu = i - ci
              call ppm(mx_d(iu-2,j,k), mx_d(iu-1,j,k), mx_d(iu,j,k), mx_d(iu+1,j,k), mx_d(iu+2,j,k), ml, dm, m6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfx_d(i,j,k) =  u_d(i,j,k) * (sum(mx_d(iu+1:i,j,k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx_d(i,j,k)
            else
              iu = i - ci + 1
              call ppm(mx_d(iu-2,j,k), mx_d(iu-1,j,k), mx_d(iu,j,k), mx_d(iu+1,j,k), mx_d(iu+2,j,k), ml, dm, m6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfx_d(i,j,k) = -u_d(i,j,k) * (sum(mx_d(i+1:iu-1,j,k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx_d(i,j,k)
            end if
          end do
        end do
      end do 
      !$omp end target teams distribute parallel do
      js = mesh%half_jds
      je = mesh%half_jde
      is = mesh%full_ids
      ie = mesh%full_ide
      !$omp target teams distribute parallel do private(i, j, k, ju, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(3)
      do k = ks, ke
        ! Along y-axis
        do j = js, je
          do i = is, ie
            if (abs(cfly_d(i,j,k)) < 1.0e-16_r8) then
              mfy_d(i,j,k) = 0
            else if (cfly_d(i,j,k) > 0) then
              ju = j
              call ppm(my_d(i,ju-2,k), my_d(i,ju-1,k), my_d(i,ju,k), my_d(i,ju+1,k), my_d(i,ju+2,k), ml, dm, m6)
              s1 = 1 - cfly_d(i,j,k)
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfy_d(i,j,k) =  v_d(i,j,k) * (ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly_d(i,j,k)
            else if (cfly_d(i,j,k) < 0) then
              ju = j + 1
              call ppm(my_d(i,ju-2,k), my_d(i,ju-1,k), my_d(i,ju,k), my_d(i,ju+1,k), my_d(i,ju+2,k), ml, dm, m6)
              s1 = 0
              s2 = -cfly_d(i,j,k)
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfy_d(i,j,k) = -v_d(i,j,k) * (ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly_d(i,j,k)
            end if
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      !$omp target update from(mfx_d, mfy_d)
      ! call mfx%DeviceToHost()
      ! call mfy%DeviceToHost()
    case ('vtx')
    end select
    end associate

    call perf_stop('hflx_ppm')
  
  end subroutine hflx_ppm

  subroutine vflx_ppm(batch, w, m, mfz)

    type(adv_batch_type     ), intent(inout) :: batch
    type(latlon_field3d_type), intent(in   ) :: w
    type(latlon_field3d_type), intent(in   ) :: m
    type(latlon_field3d_type), intent(inout) :: mfz

    integer i, j, k, ku, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6

    call perf_start('vflx_ppm')
    
    associate (mesh => m%mesh    , &
               cflz => batch%cflz)   ! in
    select case (batch%loc)
    case ('cell')
      !$omp parallel do private(i, j, k, ku, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(2)
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci - 1
              ku = min(max(ku, mesh%full_kms + 2), mesh%full_kme - 2)
              call ppm(m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2), ml, dm, m6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfz%d(i,j,k) =  w%d(i,j,k) * (sum(m%d(i,j,ku+1:k-1)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz%d(i,j,k)
            else
              ku = k - ci
              ku = min(max(ku, mesh%full_kms + 2), mesh%full_kme - 2)
              call ppm(m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2), ml, dm, m6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfz%d(i,j,k) = -w%d(i,j,k) * (sum(m%d(i,j,k:ku-1)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz%d(i,j,k)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    case ('lev')
      !$omp parallel do private(k, j, i, ci, cf, ku, ml, dm, m6, s1, s2, ds1, ds2, ds3) collapse(2)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ci = int(cflz%d(i,j,k))
            cf = cflz%d(i,j,k) - ci
            if (abs(cflz%d(i,j,k)) < 1.0e-16_r8) then
              mfz%d(i,j,k) = 0
            else if (cflz%d(i,j,k) > 0) then
              ku = k - ci
              call ppm(m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2), ml, dm, m6)
              s1 = 1 - cf
              s2 = 1
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfz%d(i,j,k) =  w%d(i,j,k) * (sum(m%d(i,j,ku+1:k)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz%d(i,j,k)
            else
              ku = k - ci + 1
              call ppm(m%d(i,j,ku-2), m%d(i,j,ku-1), m%d(i,j,ku), m%d(i,j,ku+1), m%d(i,j,ku+2), ml, dm, m6)
              s1 = 0
              s2 = -cf
              ds1 = s2    - s1
              ds2 = s2**2 - s1**2
              ds3 = s2**3 - s1**3
              mfz%d(i,j,k) = -w%d(i,j,k) * (sum(m%d(i,j,k+1:ku-1)) + ml * ds1 + 0.5_r8 * dm * ds2 + m6 * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflz%d(i,j,k)
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end select
    end associate

    call perf_stop('vflx_ppm')

  end subroutine vflx_ppm

end module ffsl_mod
