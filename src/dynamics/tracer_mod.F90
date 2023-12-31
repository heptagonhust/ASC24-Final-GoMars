! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module tracer_mod

  use flogger
  use string
  use const_mod, only: r8
  use namelist_mod
  use block_mod
  use latlon_parallel_mod
  use process_mod, only: proc
  use tracer_types_mod
  use interp_mod

  implicit none

  private

  public tracer_init_stage1
  public tracer_init_stage2
  public tracer_final
  public tracer_add
  public tracer_get_idx
  public tracer_calc_qm
  public tracer_fill_negative_values
  public ntracers
  public ntracers_water
  public nbatches
  public idx_qv
  public idx_qc, idx_nc
  public idx_qi, idx_ni
  public idx_qr, idx_nr
  public idx_qs, idx_ns
  public idx_qg
  public idx_qh
  public idx_qo3
  public idx_qso2
  public batch_names
  public batch_dts
  public tracer_batches
  public tracer_names
  public tracer_long_names
  public tracer_units
  public tracer_types
  public tracers_type
  public tracers

contains

  subroutine tracer_init_stage1()

    integer iblk

    call tracer_final()

    allocate(batch_names      (10 )); batch_names       = 'N/A'
    allocate(batch_dts        (10 )); batch_dts         = 0
    allocate(tracer_batches   (100)); tracer_batches    = 'N/A'
    allocate(tracer_names     (100)); tracer_names      = 'N/A'
    allocate(tracer_long_names(100)); tracer_long_names = 'N/A'
    allocate(tracer_units     (100)); tracer_units      = 'kg kg-1'
    allocate(tracer_types     (100)); tracer_types      = 0

    allocate(tracers(size(blocks)))
    do iblk = 1, size(blocks)
      call tracers(iblk)%init_stage1(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

  end subroutine tracer_init_stage1

  subroutine tracer_init_stage2()

    integer iblk, i

    if (ntracers == 0) return

    ! Allocate tracer arrays for each block.
    do iblk = 1, size(blocks)
      call tracers(iblk)%init_stage2(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

    do iblk = 1, size(blocks)
      ! Allocate physics tendency arrays for dynamics.
      call blocks(iblk)%aux%init_phys(blocks(iblk)%filter_mesh, blocks(iblk)%filter_halo, blocks(iblk)%mesh, blocks(iblk)%halo)
    end do

  end subroutine tracer_init_stage2

  subroutine tracer_final()

    nbatches       = 0
    ntracers       = 0
    ntracers_water = 0
    idx_qv         = 0
    idx_qc         = 0
    idx_nc         = 0
    idx_qi         = 0
    idx_ni         = 0
    idx_qr         = 0
    idx_nr         = 0
    idx_qs         = 0
    idx_ns         = 0
    idx_qg         = 0
    idx_qh         = 0
    idx_qo3        = 0

    if (allocated(batch_names      )) deallocate(batch_names      )
    if (allocated(batch_dts        )) deallocate(batch_dts        )
    if (allocated(tracer_batches   )) deallocate(tracer_batches   )
    if (allocated(tracer_names     )) deallocate(tracer_names     )
    if (allocated(tracer_long_names)) deallocate(tracer_long_names)
    if (allocated(tracer_units     )) deallocate(tracer_units     )
    if (allocated(tracer_types     )) deallocate(tracer_types     )
    if (allocated(tracers          )) deallocate(tracers          )

  end subroutine tracer_final

  subroutine tracer_add(batch_name, dt, name, long_name, units, type)

    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in), optional :: units
    integer, intent(in), optional :: type

    integer i
    logical found

    found = .false.
    do i = 1, nbatches
      if (batch_name == batch_names(i)) then
        found = .true.
        exit
      end if
    end do
    if (.not. found) then
      nbatches = nbatches + 1
      batch_names(nbatches) = batch_name
      batch_dts(nbatches) = dt
    end if

    ntracers = ntracers + 1
    tracer_batches(ntracers) = batch_name
    tracer_names(ntracers) = name
    tracer_long_names(ntracers) = long_name
    if (present(units)) tracer_units(ntracers) = units
    if (present(type)) tracer_types(ntracers) = type

    ! Set tracer indices.
    select case (name)
    case ('qv', 'Q')
      idx_qv    = ntracers; ntracers_water = ntracers_water + 1
    case ('qc', 'CLDLIQ')
      idx_qc    = ntracers; ntracers_water = ntracers_water + 1
    case ('nc', 'NUMLIQ')
      idx_nc    = ntracers
    case ('qi', 'CLDICE')
      idx_qi    = ntracers; ntracers_water = ntracers_water + 1
    case ('ni', 'NUMICE')
      idx_ni    = ntracers
    case ('qr', 'RAINQM')
      idx_qr    = ntracers; ntracers_water = ntracers_water + 1
    case ('qs', 'SNOWQM')
      idx_qs    = ntracers; ntracers_water = ntracers_water + 1
    case ('qg')
      idx_qg    = ntracers; ntracers_water = ntracers_water + 1
    case ('qh')
      idx_qh    = ntracers; ntracers_water = ntracers_water + 1
    case ('qo3')
      idx_qo3   = ntracers
    case ('qso2', 'SO2')
      idx_qso2  = ntracers
    end select

  end subroutine tracer_add

  pure integer function tracer_get_idx(name) result(res)

    character(*), intent(in) :: name

    integer i

    do i = 1, ntracers
      if (name == tracer_names(i)) then
        res = i
        return
      end if
    end do
    res = 0

  end function tracer_get_idx

  subroutine tracer_calc_qm(block)

    type(block_type), intent(in) :: block

    integer i, j, k

    if (.not. allocated(tracers)) return
    if (.not. tracers(block%id)%qm%initialized) return

    associate (mesh   => block%mesh              , &
               q      => tracers(block%id)%q     , & ! in
               qm     => tracers(block%id)%qm    , & ! out
               qm_lev => tracers(block%id)%qm_lev)   ! out
    if (idx_qv > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = q%d(i,j,k,idx_qv)
          end do
        end do
      end do
    end if
    if (idx_qc > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qc)
          end do
        end do
      end do
    end if
    if (idx_qi > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qi)
          end do
        end do
      end do
    end if
    if (idx_qr > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qr)
          end do
        end do
      end do
    end if
    if (idx_qs > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qs)
          end do
        end do
      end do
    end if
    if (idx_qg > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qg)
          end do
        end do
      end do
    end if
    if (idx_qh > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm%d(i,j,k) = qm%d(i,j,k) + q%d(i,j,k,idx_qh)
          end do
        end do
      end do
    end if
    call fill_halo(qm, west_halo=.false., south_halo=.false.)
    if (nonhydrostatic) call interp_run(qm, qm_lev)
    end associate

  end subroutine tracer_calc_qm

  subroutine tracer_fill_negative_values(block, itime, q)

    type(block_type), intent(in) :: block
    integer, intent(in) :: itime
    real(r8), intent(inout) :: q(block%filter_mesh%full_ims:block%filter_mesh%full_ime, &
                                 block%filter_mesh%full_jms:block%filter_mesh%full_jme, &
                                 block%filter_mesh%full_kms:block%filter_mesh%full_kme)

    integer i, j, k
    real(r8) neg_qm, pos_qm

    if (advection) return

    associate (mesh => block%filter_mesh, &
               dmg  => block%dstate(itime)%dmg)
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        neg_qm = 0
        pos_qm = 0
        do k = mesh%full_kds, mesh%full_kde
          if (q(i,j,k) < 0) then
            neg_qm = neg_qm + dmg%d(i,j,k) * q(i,j,k)
          else
            pos_qm = pos_qm + dmg%d(i,j,k) * q(i,j,k)
          end if
        end do
        if (neg_qm < 0) then
          if (pos_qm >= -neg_qm) then
            do k = mesh%full_kds, mesh%full_kde
              if (q(i,j,k) < 0) then
                q(i,j,k) = 0
              else
                q(i,j,k) = q(i,j,k) * (1 + neg_qm / pos_qm)
              end if
            end do
          else
            call log_warning('Negative tracer mass is larger than positive one at i=' // &
                             to_str(i) // ', j=' // to_str(j) // '!', __FILE__, __LINE__)
            do k = mesh%full_kds, mesh%full_kde
              if (q(i,j,k) < 0) then
                q(i,j,k) = 0
              end if
            end do
          end if
        end if
      end do
    end do
    end associate

  end subroutine tracer_fill_negative_values

end module tracer_mod
