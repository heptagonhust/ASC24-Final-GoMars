module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use math_mod
  use time_mod
  use block_mod
  use parallel_mod
  use tracer_mod
  use adv_batch_mod
  use ffsl_mod
  use upwind_mod
  use weno_mod
  use tvd_mod

  implicit none

  private

  public adv_init
  public adv_prepare
  public adv_run
  public adv_final
  public adv_accum_wind
  public adv_calc_mass_hflx
  public adv_calc_mass_vflx
  public adv_calc_tracer_hflx
  public adv_calc_tracer_vflx
  public adv_calc_tracer_vflx_lev
  public adv_batch_type

  public upwind1
  public upwind3
  public weno3
  public weno5
  public tvd

  interface
    subroutine calc_hflx_interface(block, batch, m, mfx, mfy, dt)
      import block_type, adv_batch_type, r8
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
    end subroutine calc_hflx_interface
    subroutine calc_vflx_interface(block, batch, m, mfz, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%full_kms:block%filter_mesh%full_kme)
      real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%half_kms:block%mesh%half_kme)
      real(r8), intent(in), optional :: dt
    end subroutine calc_vflx_interface
    subroutine calc_vflx_lev_interface(block, batch, m, mfz, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                   block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                   block%filter_mesh%half_kms:block%filter_mesh%half_kme)
      real(r8), intent(out) :: mfz(block%mesh%full_ims:block%mesh%full_ime, &
                                   block%mesh%full_jms:block%mesh%full_jme, &
                                   block%mesh%full_kms:block%mesh%full_kme)
      real(r8), intent(in), optional :: dt
    end subroutine calc_vflx_lev_interface
  end interface

  procedure(calc_hflx_interface    ), pointer :: adv_calc_mass_hflx       => null()
  procedure(calc_vflx_interface    ), pointer :: adv_calc_mass_vflx       => null()
  procedure(calc_hflx_interface    ), pointer :: adv_calc_tracer_hflx     => null()
  procedure(calc_vflx_interface    ), pointer :: adv_calc_tracer_vflx     => null()
  procedure(calc_vflx_lev_interface), pointer :: adv_calc_tracer_vflx_lev => null()

contains

  subroutine adv_init()

    integer iblk, ibat, itra, n, idx(1000)

    call adv_final()

    select case (adv_scheme)
    case ('ffsl')
      call ffsl_init()
      adv_calc_mass_hflx       => ffsl_calc_mass_hflx
      adv_calc_mass_vflx       => ffsl_calc_mass_vflx
      adv_calc_tracer_hflx     => ffsl_calc_tracer_hflx
      adv_calc_tracer_vflx     => ffsl_calc_tracer_vflx
      adv_calc_tracer_vflx_lev => ffsl_calc_tracer_vflx_lev
    case default
      call log_error('Invalid adv_scheme ' // trim(adv_scheme) // '!', pid=proc%id)
    end select

    call tvd_init()

    ! Initialize advection batches.
    do iblk = 1, size(blocks)
      if (.not. advection) then
        call blocks(iblk)%adv_batch_pt%init(           &
          blocks(iblk)%filter_mesh, blocks(iblk)%mesh, &
          'cell', 'pt', dt_dyn, dynamic=.true.)
      end if
    end do

    if (nbatches == 0) then
      call log_warning('No advection batches have been defined yet!', pid=proc%id)
      return
    end if

    do iblk = 1, size(blocks)
      allocate(blocks(iblk)%adv_batches(nbatches))
      do ibat = 1, nbatches
        n = 0
        do itra = 1, ntracers
          if (batch_names(ibat) == tracer_batches(itra)) then
            n = n + 1
            idx(n) = itra
          end if
        end do
        call blocks(iblk)%adv_batches(ibat)%init(    &
          blocks(iblk)%filter_mesh, blocks(iblk)%mesh, &
          'cell', batch_names(ibat), batch_dts(ibat), dynamic=.false., idx=idx(1:n))
      end do
    end do

    if (proc%is_root()) then
      call log_notice('There are ' // to_str(size(blocks(1)%adv_batches)) // ' advection batches.')
      do ibat = 1, size(blocks(1)%adv_batches)
        write(*, *) '- ', trim(blocks(1)%adv_batches(ibat)%name), int(blocks(1)%adv_batches(ibat)%dt)
      end do
    end if

  end subroutine adv_init

  subroutine adv_prepare(itime)

    integer, intent(in) :: itime

    integer iblk, m

    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)              , &
                 dstate => blocks(iblk)%dstate(itime))
      if (allocated(block%adv_batches)) then
        do m = 1, size(block%adv_batches)
          call block%adv_batches(m)%copy_old_m(dstate%dmg)
        end do
      end if
      call block%adv_batch_pt%copy_old_m(dstate%dmg)
      end associate
    end do
    if (.not. restart) call adv_accum_wind(itime)

  end subroutine adv_prepare

  subroutine adv_run(itime)

    integer, intent(in) :: itime

    integer iblk, i, j, k, l, m
    real(r8), pointer :: q_new(:,:,:)
    real(r8), allocatable :: q_old(:,:,:)
    real(r8), allocatable :: work(:,:), pole(:)
    real(r8) qm0, qm1, qm2, qm0_half
    real(r8), dimension(nlev) :: a, b, c, r

    call adv_accum_wind(itime)

    if (.not. allocated(blocks(1)%adv_batches)) return

    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)                  , &
                 dstate => blocks(iblk)%dstate(itime)    , &
                 mesh   => blocks(iblk)%filter_mesh      , &
                 m_new  => blocks(iblk)%dstate(itime)%dmg)
      allocate(q_old(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      allocate(work(mesh%full_ids:mesh%full_ide,mesh%full_nlev))
      allocate(pole(mesh%full_nlev))
      do m = 1, size(block%adv_batches)
        if (time_is_alerted(block%adv_batches(m)%name)) then
          associate (batch => block%adv_batches(m), &
                     is    => mesh%full_ims       , &
                     ie    => mesh%full_ime       , &
                     js    => mesh%full_jms       , &
                     je    => mesh%full_jme       , &
                     ks    => mesh%full_kms       , &
                     ke    => mesh%full_kme       )
          do l = 1, block%adv_batches(m)%ntracers
            q_old(is:ie,js:je,ks:ke) =  tracers(block%id)%q(:,:,:,batch%idx(l))
            q_new(is:ie,js:je,ks:ke) => tracers(block%id)%q(:,:,:,batch%idx(l))
            associate (m_old   => batch%old_m  , & ! inout
                       we_imp  => batch%we_imp , & ! in
                       qmf_lon => batch%qmf_lon, & ! working array
                       qmf_lat => batch%qmf_lat, & ! working array
                       qmf_lev => batch%qmf_lev)   ! working array
            ! Calculate tracer mass flux.
            call adv_calc_tracer_hflx(block, batch, q_old, qmf_lon, qmf_lat)
            call fill_halo(block%halo, qmf_lon, full_lon=.false., full_lat=.true., full_lev=.true., &
                           south_halo=.false., north_halo=.false., east_halo=.false.)
            call fill_halo(block%halo, qmf_lat, full_lon=.true., full_lat=.false., full_lev=.true., &
                           north_halo=.false.,  west_halo=.false., east_halo=.false.)
            ! Update tracer mixing ratio.
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
                do i = mesh%full_ids, mesh%full_ide
                  q_new(i,j,k) = m_old(i,j,k) * q_old(i,j,k) - ( &
                    (                                            &
                      qmf_lon(i  ,j,k) -                         &
                      qmf_lon(i-1,j,k)                           &
                    ) * mesh%le_lon(j) + (                       &
                      qmf_lat(i,j  ,k) * mesh%le_lat(j  ) -      &
                      qmf_lat(i,j-1,k) * mesh%le_lat(j-1)        &
                    )                                            &
                  ) / mesh%area_cell(j) * dt_adv
                end do
              end do
            end do
            if (mesh%has_south_pole()) then
              j = mesh%full_jds
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  work(i,k) = qmf_lat(i,j,k)
                end do
              end do
              call zonal_sum(proc%zonal_circle, work, pole)
              pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j) * dt_adv
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  q_new(i,j,k) = m_old(i,j,k) * q_old(i,j,k) - pole(k)
                end do
              end do
            end if
            if (mesh%has_north_pole()) then
              j = mesh%full_jde
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  work(i,k) = qmf_lat(i,j-1,k)
                end do
              end do
              call zonal_sum(proc%zonal_circle, work, pole)
              pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j) * dt_adv
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  q_new(i,j,k) = m_old(i,j,k) * q_old(i,j,k) + pole(k)
                end do
              end do
            end if
            ! Fill possible negative values.
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  if (q_new(i,j,k) < 0) then
                    qm0 = q_new(i,j,k)
                    qm1 = merge(q_new(i,j,k-1), 0.0_r8, k > mesh%full_kds)
                    qm2 = merge(q_new(i,j,k+1), 0.0_r8, k < mesh%full_kde)
                    qm0_half = 0.5_r8 * qm0
                    if (qm1 >= qm0_half .and. qm2 >= qm0_half) then
                      if (qm1 > 0) q_new(i,j,k-1) = qm1 - qm0_half
                      if (qm2 > 0) q_new(i,j,k+1) = qm2 - qm0_half
                    else if (qm1 > qm0) then
                      if (qm1 > 0) q_new(i,j,k-1) = qm1 - qm0
                    else if (qm2 > qm0) then
                      if (qm2 > 0) q_new(i,j,k+1) = qm2 - qm0
                    else
                      call log_error('Negative tracer!')
                    end if
                    q_new(i,j,k) = 0
                  end if
                end do
              end do
            end do
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  q_new(i,j,k) = q_new(i,j,k) / m_new(i,j,k)
                end do
              end do
            end do
            ! Set upper and lower boundary conditions.
            do k = mesh%full_kms, mesh%full_kds - 1
              q_new(:,:,k) = q_new(:,:,mesh%full_kds)
            end do
            do k = mesh%full_kde + 1, mesh%full_kme
              q_new(:,:,k) = q_new(:,:,mesh%full_kde)
            end do
            call adv_calc_tracer_vflx(block, block%adv_batches(m), q_new, qmf_lev)
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  q_new(i,j,k) = q_new(i,j,k) * m_new(i,j,k) - (qmf_lev(i,j,k+1) - qmf_lev(i,j,k)) * dt_adv
                end do
              end do
            end do
            if (use_ieva) then
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  a(1) = 0
                  c(mesh%full_kde) = 0
                  do k = mesh%full_kds + 1, mesh%full_kde
                    a(k)   = -dt_adv * we_imp(i,j,k) * (1 + sign(1.0_r8, we_imp(i,j,k)))
                    c(k-1) =  dt_adv * we_imp(i,j,k) * (1 - sign(1.0_r8, we_imp(i,j,k)))
                  end do
                  do k = mesh%full_kds, mesh%full_kde
                    b(k) = dt_adv * we_imp(i,j,k+1) * (1 + sign(1.0_r8, we_imp(i,j,k+1))) - &
                           dt_adv * we_imp(i,j,k  ) * (1 - sign(1.0_r8, we_imp(i,j,k  ))) + &
                           2 * m_new(i,j,k)
                    r(k) = 2 * q_new(i,j,k)
                  end do
                  call triiag_thomas(a, b, c, r, q_new(i,j,mesh%full_kds:mesh%full_kde))
                end do
              end do
            end if
            if (.not. use_ieva) then
              do k = mesh%full_kds, mesh%full_kde
                do j = mesh%full_jds, mesh%full_jde
                  do i = mesh%full_ids, mesh%full_ide
                    q_new(i,j,k) = q_new(i,j,k) / m_new(i,j,k)
                  end do
                end do
              end do
            end if
            call fill_halo(block%filter_halo, q_new, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
            end associate
          end do
          end associate
          if (block%adv_batches(m)%name == 'moist') call tracer_calc_qm(block)
        end if
        call block%adv_batches(m)%copy_old_m(m_new)
      end do
      deallocate(q_old, work, pole)
      end associate
    end do

  end subroutine adv_run

  subroutine adv_accum_wind(itime)

    integer, intent(in) :: itime

    integer iblk, l

    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                      , &
                 u_lon   => blocks(iblk)%dstate(itime)%u_lon  , &
                 v_lat   => blocks(iblk)%dstate(itime)%v_lat  , &
                 mfx_lon => blocks(iblk)%dstate(itime)%mfx_lon, &
                 mfy_lat => blocks(iblk)%dstate(itime)%mfy_lat, &
                 we_lev  => blocks(iblk)%dstate(itime)%we_lev , &
                 dmg_lev => blocks(iblk)%dstate(itime)%dmg_lev)
      if (allocated(block%adv_batches)) then
        do l = 1, size(block%adv_batches)
          select case (block%adv_batches(l)%loc)
          case ('cell')
            call block%adv_batches(l)%accum_uv_cell(u_lon, v_lat)
            call block%adv_batches(l)%accum_mf_cell(mfx_lon, mfy_lat)
            if (global_mesh%full_nlev > 1) then
              call block%adv_batches(l)%accum_we_lev(we_lev, dmg_lev)
            end if
          end select
        end do
      end if
      end associate
    end do

  end subroutine adv_accum_wind

  subroutine adv_final()

  end subroutine adv_final

end module adv_mod
