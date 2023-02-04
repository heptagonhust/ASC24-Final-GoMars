module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use parallel_mod
  use adv_batch_mod
  use ffsl_mod
  use upwind_mod
  use weno_mod
  use tvd_mod
  use filter_mod

  implicit none

  private

  public adv_init
  public adv_prepare
  public adv_run
  public adv_final
  public adv_add_tracer
  public adv_allocate_tracers
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

  interface adv_allocate_tracers
    module procedure adv_allocate_tracers_1
    module procedure adv_allocate_tracers_2
  end interface adv_allocate_tracers

  procedure(calc_hflx_interface    ), pointer :: adv_calc_mass_hflx       => null()
  procedure(calc_vflx_interface    ), pointer :: adv_calc_mass_vflx       => null()
  procedure(calc_hflx_interface    ), pointer :: adv_calc_tracer_hflx     => null()
  procedure(calc_vflx_interface    ), pointer :: adv_calc_tracer_vflx     => null()
  procedure(calc_vflx_lev_interface), pointer :: adv_calc_tracer_vflx_lev => null()

  integer ntracer
  character(30), allocatable :: batch_names(:)
  character(30), allocatable :: tracer_names(:)
  character(30), allocatable :: tracer_long_names(:)
  character(30), allocatable :: tracer_units(:)
  real(r8), allocatable :: tracer_dt(:)

contains

  subroutine adv_init()

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

    ntracer = 0
    allocate(batch_names      (1000))
    allocate(tracer_names     (1000))
    allocate(tracer_long_names(1000))
    allocate(tracer_units     (1000))
    allocate(tracer_dt        (1000))

    call tvd_init()

  end subroutine adv_init

  subroutine adv_prepare(itime)

    integer, intent(in) :: itime

    integer iblk, m

    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)              , &
                 dstate => blocks(iblk)%dstate(itime))
      if (allocated(block%adv_batches)) then
        do m = 1, size(block%adv_batches)
          call block%adv_batches(m)%copy_old_m(dstate%m)
        end do
      end if
      call block%adv_batch_pt%copy_old_m(dstate%m)
      if (.not. restart) call adv_accum_wind(block, itime)
      end associate
    end do

  end subroutine adv_prepare

  subroutine adv_run(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer i, j, k, l, m
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev), qm0, qm1, qm2, qm0_half

    call adv_accum_wind(block, itime)

    if (.not. allocated(block%adv_batches)) return

    do m = 1, size(block%adv_batches)
      if (time_is_alerted(block%adv_batches(m)%name)) then
        do l = 1, size(block%adv_batches(m)%tracer_names)
          associate (mesh    => block%mesh                  , &
                     old     => block%adv_batches(m)%old    , &
                     new     => block%adv_batches(m)%new    , &
                     old_m   => block%adv_batches(m)%old_m  , & ! inout
                     q       => block%adv_batches(m)%q      , & ! inout
                     qmf_lon => block%adv_batches(m)%qmf_lon, & ! working array
                     qmf_lat => block%adv_batches(m)%qmf_lat, & ! working array
                     qmf_lev => block%adv_batches(m)%qmf_lev)   ! working array
          ! Calculate tracer mass flux.
          call adv_calc_tracer_hflx(block, block%adv_batches(m), q(:,:,:,l,old), qmf_lon, qmf_lat)
          call fill_halo(block%halo, qmf_lon, full_lon=.false., full_lat=.true., full_lev=.true., &
                         south_halo=.false., north_halo=.false., east_halo=.false.)
          call fill_halo(block%halo, qmf_lat, full_lon=.true., full_lat=.false., full_lev=.true., &
                         north_halo=.false.,  west_halo=.false., east_halo=.false.)
          ! Update tracer mixing ratio.
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
              do i = mesh%full_ids, mesh%full_ide
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) - ( &
                  (                                                &
                    qmf_lon(i  ,j,k) -                             &
                    qmf_lon(i-1,j,k)                               &
                  ) * mesh%le_lon(j) + (                           &
                    qmf_lat(i,j  ,k) * mesh%le_lat(j  ) -          &
                    qmf_lat(i,j-1,k) * mesh%le_lat(j-1)            &
                  )                                                &
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
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) - pole(k)
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
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) + pole(k)
              end do
            end do
          end if
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                q(i,j,k,l,new) = q(i,j,k,l,new) / block%dstate(itime)%m(i,j,k)
              end do
            end do
          end do
          ! Set upper and lower boundary conditions.
          do k = mesh%full_kms, mesh%full_kds - 1
            q(:,:,k,l,new) = q(:,:,mesh%full_kds,l,new)
          end do
          do k = mesh%full_kde + 1, mesh%full_kme
            q(:,:,k,l,new) = q(:,:,mesh%full_kde,l,new)
          end do
          call adv_calc_tracer_vflx(block, block%adv_batches(m), q(:,:,:,l,new), qmf_lev)
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                q(i,j,k,l,new) = q(i,j,k,l,new) * block%dstate(itime)%m(i,j,k) - (qmf_lev(i,j,k+1) - qmf_lev(i,j,k)) * dt_adv
              end do
            end do
          end do
          ! Fill possible negative values.
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                if (q(i,j,k,l,new) < 0) then
                  qm0 = q(i,j,k  ,l,new)
                  qm1 = merge(q(i,j,k-1,l,new), 0.0_r8, k > mesh%full_kds)
                  qm2 = merge(q(i,j,k+1,l,new), 0.0_r8, k < mesh%full_kde)
                  qm0_half = 0.5_r8 * qm0
                  if (qm1 >= qm0_half .and. qm2 >= qm0_half) then
                    if (qm1 > 0) q(i,j,k-1,l,new) = qm1 - qm0_half
                    if (qm2 > 0) q(i,j,k+1,l,new) = qm2 - qm0_half
                  else if (qm1 > qm0) then
                    if (qm1 > 0) q(i,j,k-1,l,new) = qm1 - qm0
                  else if (qm2 > qm0) then
                    if (qm2 > 0) q(i,j,k+1,l,new) = qm2 - qm0
                  else
                    call log_error('Negative tracer!')
                  end if
                  q(i,j,k,l,new) = 0
                end if
              end do
            end do
          end do
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                q(i,j,k,l,new) = q(i,j,k,l,new) / block%dstate(itime)%m(i,j,k)
              end do
            end do
          end do
          call fill_halo(block%filter_halo, q(:,:,:,l,new), full_lon=.true., full_lat=.true., full_lev=.true.)
          end associate
        end do
        i = block%adv_batches(m)%old
        block%adv_batches(m)%old = block%adv_batches(m)%new
        block%adv_batches(m)%new = i
      end if
      call block%adv_batches(m)%copy_old_m(block%dstate(itime)%m)
    end do

  end subroutine adv_run

  subroutine adv_add_tracer(batch_name, dt, name, long_name, units)

    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    character(*), intent(in) :: name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units

    ntracer = ntracer + 1
    if (ntracer > size(batch_names)) then
      call log_error('Insufficient character array!', __FILE__, __LINE__, pid=proc%id)
    end if
    batch_names (ntracer) = batch_name
    tracer_names(ntracer) = name
    tracer_dt   (ntracer) = dt
    if (present(long_name)) tracer_long_names(ntracer) = long_name
    if (present(units    )) tracer_units     (ntracer) = units

  end subroutine adv_add_tracer

  subroutine adv_allocate_tracers_1(block)

    type(block_type), intent(inout) :: block

    integer nbatch, nbatch_tracer, i, j, k
    logical found
    character(30) unique_batch_names(100)
    real(r8) unique_tracer_dt(100)

    if (.not. advection) then
      call block%adv_batch_pt%init(block%filter_mesh, block%mesh, 'cell', 'pt', dt_dyn, dynamic=.true.)
    end if

    nbatch = 0
    do i = 1, ntracer
      found = .false.
      do j = 1, i - 1
        if (i /= j .and. batch_names(i) == batch_names(j)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        nbatch = nbatch + 1
        unique_batch_names(nbatch) = batch_names(i)
        unique_tracer_dt  (nbatch) = tracer_dt  (i)
      end if
    end do
    if (nbatch == 0) return
    if (allocated(block%adv_batches)) then
      call log_error('Advection batches have already been alllocated!', pid=proc%id)
    end if
    if (proc%is_root()) then
      call log_notice('There are ' // to_str(nbatch) // ' advection batches.')
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_tracer_dt(i))
      end do
    end if

    ! Initialize advection batches in block objects and allocate tracer arrays in dstate objects.
    allocate(block%adv_batches(nbatch))
    do i = 1, nbatch
      call block%adv_batches(i)%init(block%filter_mesh, block%mesh, 'cell', unique_batch_names(i), unique_tracer_dt(i), dynamic=.false.)
    end do

    ! Record tracer information in advection batches.
    do i = 1, nbatch
      nbatch_tracer = 0
      do j = 1, ntracer
        if (batch_names(j) == block%adv_batches(i)%name) then
          nbatch_tracer = nbatch_tracer + 1
        end if
      end do
      call block%adv_batches(i)%allocate_tracers(nbatch_tracer)
      associate (filter_mesh => block%filter_mesh, mesh => block%mesh)
      allocate(block%adv_batches(i)%q(filter_mesh%full_ims:filter_mesh%full_ime, &
                                      filter_mesh%full_jms:filter_mesh%full_jme, &
                                      filter_mesh%full_kms:filter_mesh%full_kme,nbatch_tracer,2))
      allocate(block%adv_batches(i)%qmf_lon(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mesh%full_kms:mesh%full_kme))
      allocate(block%adv_batches(i)%qmf_lat(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mesh%full_kms:mesh%full_kme))
      allocate(block%adv_batches(i)%qmf_lev(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mesh%half_kms:mesh%half_kme))
      end associate
      k = 0
      do j = 1, ntracer
        if (batch_names(j) == block%adv_batches(i)%name) then
          k = k + 1
          block%adv_batches(i)%tracer_names     (k) = tracer_names     (j)
          block%adv_batches(i)%tracer_long_names(k) = tracer_long_names(j)
          block%adv_batches(i)%tracer_units     (k) = tracer_units     (j)
        end if
      end do
    end do

  end subroutine adv_allocate_tracers_1

  subroutine adv_allocate_tracers_2(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call adv_allocate_tracers_1(blocks(iblk))
    end do

  end subroutine adv_allocate_tracers_2

  subroutine adv_accum_wind(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer l

    if (allocated(block%adv_batches)) then
      do l = 1, size(block%adv_batches)
        select case (block%adv_batches(l)%loc)
        case ('cell')
          call block%adv_batches(l)%accum_uv_cell( &
            block%dstate(itime)%u_lon            , &
            block%dstate(itime)%v_lat            )
          call block%adv_batches(l)%accum_mf_cell( &
            block%dstate(itime)%mfx_lon          , &
            block%dstate(itime)%mfy_lat          )
          if (global_mesh%full_nlev > 1) then
            call block%adv_batches(l)%accum_we_lev( &
              block%dstate(itime)%we_lev          , &
              block%dstate(itime)%m_lev           )
          end if
        end select
      end do
    end if

  end subroutine adv_accum_wind

  subroutine adv_final()

    ntracer = 0
    if (allocated(batch_names      )) deallocate(batch_names      )
    if (allocated(tracer_names     )) deallocate(tracer_names     )
    if (allocated(tracer_long_names)) deallocate(tracer_long_names)
    if (allocated(tracer_units     )) deallocate(tracer_units     )
    if (allocated(tracer_dt        )) deallocate(tracer_dt        )

  end subroutine adv_final

end module adv_mod
