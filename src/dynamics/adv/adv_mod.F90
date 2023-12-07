! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use math_mod
  use time_mod
  use block_mod
  use latlon_field_types_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use tracer_mod
  use adv_batch_mod
  use ffsl_mod
  use upwind_mod
  use weno_mod
  use tvd_mod
  use physics_mod

  implicit none

  private

  public adv_init
  public adv_prepare
  public adv_run
  public adv_final
  public adv_fill_vhalo
  public adv_accum_wind
  public adv_calc_mass_hflx
  public adv_calc_mass_vflx
  public adv_calc_tracer_hflx
  public adv_calc_tracer_vflx
  public adv_batch_type

  public upwind1
  public upwind3
  public upwind5
  public weno3
  public weno5
  public tvd

  interface
    subroutine calc_hflx_interface(batch, m, mfx, mfy, dt)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfx
      type(latlon_field3d_type), intent(inout) :: mfy
      real(r8), intent(in), optional :: dt
    end subroutine calc_hflx_interface
    subroutine calc_vflx_interface(batch, m, mfz, dt)
      import adv_batch_type, latlon_field3d_type, r8
      type(adv_batch_type     ), intent(inout) :: batch
      type(latlon_field3d_type), intent(in   ) :: m
      type(latlon_field3d_type), intent(inout) :: mfz
      real(r8), intent(in), optional :: dt
    end subroutine calc_vflx_interface
  end interface

  procedure(calc_hflx_interface), pointer :: adv_calc_mass_hflx   => null()
  procedure(calc_vflx_interface), pointer :: adv_calc_mass_vflx   => null()
  procedure(calc_hflx_interface), pointer :: adv_calc_tracer_hflx => null()
  procedure(calc_vflx_interface), pointer :: adv_calc_tracer_vflx => null()

contains

  subroutine adv_init()

    integer iblk, ibat, itra, n, idx(1000)

    call adv_final()

    select case (adv_scheme)
    case ('ffsl')
      call ffsl_init()
      adv_calc_mass_hflx   => ffsl_calc_mass_hflx
      adv_calc_mass_vflx   => ffsl_calc_mass_vflx
      adv_calc_tracer_hflx => ffsl_calc_tracer_hflx
      adv_calc_tracer_vflx => ffsl_calc_tracer_vflx
    case default
      call log_error('Invalid adv_scheme ' // trim(adv_scheme) // '!', pid=proc%id)
    end select

    call tvd_init()

    ! Initialize advection batches.
    do iblk = 1, size(blocks)
      if (.not. advection) then
        call blocks(iblk)%adv_batch_pt%init(                  &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
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
        call blocks(iblk)%adv_batches(ibat)%init(             &
          blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, &
          blocks(iblk)%mesh, blocks(iblk)%halo              , &
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

    if (.not. restart) then
      do iblk = 1, size(blocks)
        associate (dmg => blocks(iblk)%dstate(itime)%dmg)
        if (allocated(blocks(iblk)%adv_batches)) then
          do m = 1, size(blocks(iblk)%adv_batches)
            call blocks(iblk)%adv_batches(m)%copy_old_m(dmg)
          end do
        end if
        end associate
      end do
      call adv_accum_wind(itime)
    end if

  end subroutine adv_prepare

  subroutine adv_run(itime)

    integer, intent(in) :: itime

    integer iblk, i, j, k, l, m, idx
    type(latlon_field3d_type) q_old, q_new
    real(r8), allocatable :: work(:,:), pole(:)

    if (.not. allocated(blocks(1)%adv_batches)) return

    call adv_accum_wind(itime)

    do iblk = 1, size(blocks)
      associate (block     => blocks(iblk)                  , &
                 dstate    => blocks(iblk)%dstate(itime)    , &
                 mesh      => blocks(iblk)%filter_mesh      , &
                 m_new     => blocks(iblk)%dstate(itime)%dmg)
      allocate(work(mesh%full_ids:mesh%full_ide,mesh%full_nlev))
      allocate(pole(mesh%full_nlev))
      do m = 1, size(block%adv_batches)
        if (time_is_alerted(block%adv_batches(m)%name)) then
          if (m == 1 .and. pdc_type == 2) call physics_update_dynamics(block, itime, dt_adv)
          associate (batch => block%adv_batches(m))
          call q_old%link(batch%dmf) ! Borrow array.
          do l = 1, block%adv_batches(m)%ntracers
            idx = batch%idx(l)
            call q_new%link(tracers(iblk)%q, idx)
            q_old%d = q_new%d
            associate (m_old   => batch%old_m  , & ! inout
                       qmf_lon => batch%qmf_lon, & ! working array
                       qmf_lat => batch%qmf_lat, & ! working array
                       qmf_lev => batch%qmf_lev)   ! working array
            ! Calculate tracer mass flux.
            call adv_calc_tracer_hflx(batch, q_old, qmf_lon, qmf_lat)
            call fill_halo(qmf_lon, south_halo=.false., north_halo=.false., east_halo=.false.)
            call fill_halo(qmf_lat, north_halo=.false.,  west_halo=.false., east_halo=.false.)
            ! Update tracer mixing ratio.
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = (m_old%d(i,j,k) * q_old%d(i,j,k) - ( &
                    (                                                   &
                      qmf_lon%d(i  ,j,k) -                              &
                      qmf_lon%d(i-1,j,k)                                &
                    ) * mesh%le_lon(j) + (                              &
                      qmf_lat%d(i,j  ,k) * mesh%le_lat(j  ) -           &
                      qmf_lat%d(i,j-1,k) * mesh%le_lat(j-1)             &
                    )                                                   &
                  ) / mesh%area_cell(j) * dt_adv) / m_new%d(i,j,k)
                end do
              end do
            end do
            if (mesh%has_south_pole()) then
              j = mesh%full_jds
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  work(i,k) = qmf_lat%d(i,j,k)
                end do
              end do
              call zonal_sum(proc%zonal_circle, work, pole)
              pole = pole * mesh%le_lat(j) / global_mesh%full_nlon / mesh%area_cell(j) * dt_adv
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = (m_old%d(i,j,k) * q_old%d(i,j,k) - pole(k)) / m_new%d(i,j,k)
                end do
              end do
            end if
            if (mesh%has_north_pole()) then
              j = mesh%full_jde
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  work(i,k) = qmf_lat%d(i,j-1,k)
                end do
              end do
              call zonal_sum(proc%zonal_circle, work, pole)
              pole = pole * mesh%le_lat(j-1) / global_mesh%full_nlon / mesh%area_cell(j) * dt_adv
              do k = mesh%full_kds, mesh%full_kde
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = (m_old%d(i,j,k) * q_old%d(i,j,k) + pole(k)) / m_new%d(i,j,k)
                end do
              end do
            end if
            call adv_fill_vhalo(q_new)
            call adv_calc_tracer_vflx(block%adv_batches(m), q_new, qmf_lev)
            do k = mesh%full_kds, mesh%full_kde
              do j = mesh%full_jds, mesh%full_jde
                do i = mesh%full_ids, mesh%full_ide
                  q_new%d(i,j,k) = q_new%d(i,j,k) - (qmf_lev%d(i,j,k+1) - qmf_lev%d(i,j,k)) * dt_adv / m_new%d(i,j,k)
                end do
              end do
            end do
            if (pdc_type == 1 .or. pdc_type == 2) then
              call physics_update_tracers(block, itime, dt_adv, batch%idx(l))
            else
              call tracer_fill_negative_values(block, itime, q_new%d)
              call fill_halo(q_new)
            end if
            end associate
          end do
          end associate
          call block%adv_batches(m)%copy_old_m(m_new)
        end if
      end do
      deallocate(work, pole)
      call tracer_calc_qm(block)
      end associate
    end do

  end subroutine adv_run

  subroutine adv_fill_vhalo(f)

    type(latlon_field3d_type), intent(inout) :: f

    integer i, j, k

    ! Set upper and lower boundary conditions.
    do k = f%mesh%full_kds - 1, f%mesh%full_kms, -1
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,mesh%full_kds)
          ! f%d(i,j,k) = 2 * f%d(i,j,k+1) - f%d(i,j,k+2)
          f%d(i,j,k) = 3 * f%d(i,j,k+1) - 3 * f%d(i,j,k+2) + f%d(i,j,k+3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k+1) - 6 * f%d(i,j,k+2) + 4 * f%d(i,j,k+3) - f%d(i,j,k+4)
        end do
      end do
    end do
    do k = f%mesh%full_kde + 1, f%mesh%full_kme
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,mesh%full_kde)
          ! f%d(i,j,k) = 2 * f%d(i,j,k-1) - f%d(i,j,k-2)
          f%d(i,j,k) = 3 * f%d(i,j,k-1) - 3 * f%d(i,j,k-2) + f%d(i,j,k-3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k-1) - 6 * f%d(i,j,k-2) + 4 * f%d(i,j,k-3) - f%d(i,j,k-4)
        end do
      end do
    end do

  end subroutine adv_fill_vhalo

  subroutine adv_accum_wind(itime)

    integer, intent(in) :: itime

    integer iblk, l

    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                      , &
                 dmg_lon => blocks(iblk)%aux%dmg_lon          , &
                 dmg_lat => blocks(iblk)%aux%dmg_lat          , &
                 dmg_lev => blocks(iblk)%dstate(itime)%dmg_lev, &
                 mfx_lon => blocks(iblk)%aux%mfx_lon          , &
                 mfy_lat => blocks(iblk)%aux%mfy_lat          )
      if (allocated(block%adv_batches)) then
        do l = 1, size(block%adv_batches)
          select case (block%adv_batches(l)%loc)
          case ('cell')
            call block%adv_batches(l)%accum_wind(dmg_lon, dmg_lat, dmg_lev, mfx_lon, mfy_lat)
          end select
        end do
      end if
      end associate
    end do

  end subroutine adv_accum_wind

  subroutine adv_final()

  end subroutine adv_final

end module adv_mod
