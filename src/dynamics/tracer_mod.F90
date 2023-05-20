module tracer_mod

  use const_mod, only: r8
  use namelist_mod, only: mp_scheme, dt_adv
  use block_mod
  use parallel_mod
  use tracer_types_mod

  implicit none

  private

  public tracer_init
  public tracer_final
  public tracer_add
  public tracer_add_moist
  public tracer_allocate
  public tracer_get_idx
  public tracer_get_array
  public tracer_get_array_qm
  public tracer_calc_qm
  public ntracers
  public ntracers_water
  public nbatches
  public idx_qv
  public idx_qc
  public idx_qi
  public idx_qr
  public idx_qs
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
  public tracers

  interface tracer_get_array
    module procedure tracer_get_array_idx
    module procedure tracer_get_array_name
  end interface tracer_get_array

contains

  subroutine tracer_init()

    call tracer_final()

    allocate(batch_names      (10 )); batch_names       = 'N/A'
    allocate(batch_dts        (10 )); batch_dts         = 0
    allocate(tracer_batches   (100)); tracer_batches    = 'N/A'
    allocate(tracer_names     (100)); tracer_names      = 'N/A'
    allocate(tracer_long_names(100)); tracer_long_names = 'N/A'
    allocate(tracer_units     (100)); tracer_units      = 'kg kg-1'
    allocate(tracer_types     (100)); tracer_types      = 0

  end subroutine tracer_init

  subroutine tracer_final()

    nbatches       = 0
    ntracers       = 0
    ntracers_water = 0
    idx_qv         = 0
    idx_qc         = 0
    idx_qi         = 0
    idx_qr         = 0
    idx_qs         = 0
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

  subroutine tracer_add_moist()

      call tracer_add('moist', dt_adv, 'qv', 'water vapor' , 'kg kg-1', type=0)
    if (mp_scheme /= 'N/A') then
      call tracer_add('moist', dt_adv, 'qc', 'cloud liquid', 'kg kg-1', type=0)
      call tracer_add('moist', dt_adv, 'qi', 'cloud ice'   , 'kg kg-1', type=0)
      call tracer_add('moist', dt_adv, 'qr', 'rain'        , 'kg kg-1', type=0)
      call tracer_add('moist', dt_adv, 'qs', 'snow'        , 'kg kg-1', type=0)
      call tracer_add('moist', dt_adv, 'qg', 'graupel'     , 'kg kg-1', type=0)
      call tracer_add('moist', dt_adv, 'qh', 'hail'        , 'kg kg-1', type=0)
    end if

  end subroutine tracer_add_moist

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
    if (present(type )) tracer_types(ntracers) = type

  end subroutine tracer_add

  subroutine tracer_allocate()

    integer iblk, ibat, i

    if (nbatches == 0) return

    ! Set tracer indices.
    do i = 1, ntracers
      select case (tracer_names(i))
      case ('qv', 'Q')
        idx_qv    = i; ntracers_water = ntracers_water + 1
      case ('qc', 'CLDLIQ')
        idx_qc    = i; ntracers_water = ntracers_water + 1
      case ('qi', 'CLDICE')
        idx_qi    = i; ntracers_water = ntracers_water + 1
      case ('qr', 'RAINQM')
        idx_qr    = i; ntracers_water = ntracers_water + 1
      case ('qs', 'SNOWQM')
        idx_qs    = i; ntracers_water = ntracers_water + 1
      case ('qg')
        idx_qg    = i; ntracers_water = ntracers_water + 1
      case ('qh')
        idx_qh    = i; ntracers_water = ntracers_water + 1
      case ('qo3')
        idx_qo3   = i
      case ('qso2', 'SO2')
        idx_qso2  = i
      end select
    end do

    ! Allocate tracer arrays for each block.
    allocate(tracers(size(blocks)))
    do iblk = 1, size(blocks)
      call tracers(iblk)%init(blocks(iblk)%mesh, blocks(iblk)%filter_mesh)
    end do

    do iblk = 1, size(blocks)
      ! Allocate tracer arrays in physics state and tendency.
      call blocks(iblk)%pstate%init(blocks(iblk)%mesh)
      call blocks(iblk)%ptend%init(blocks(iblk)%mesh)
      ! Allocate tendency arrays in dynamics tendency.
      do i = 1, size(blocks(iblk)%dtend)
        call blocks(iblk)%dtend(i)%init_phys()
      end do
    end do

  end subroutine tracer_allocate

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

  subroutine tracer_get_array_idx(iblk, idx, q, file, line)

    integer, intent(in) :: iblk
    integer, intent(in) :: idx
    real(r8), intent(out), pointer :: q(:,:,:)
    character(*), intent(in) :: file
    integer, intent(in) :: line

    if (idx < 1) then
      call log_error('Failed to get tracer array!', file, line, pid=proc%id)
    end if
    associate (mesh => tracers(iblk)%filter_mesh)
    ! NOTE: q is on filter_mesh.
    q(mesh%full_ims:mesh%full_ime, &
      mesh%full_jms:mesh%full_jme, &
      mesh%full_kms:mesh%full_kme) => tracers(iblk)%q(:,:,:,idx)
    end associate

  end subroutine tracer_get_array_idx

  subroutine tracer_get_array_name(iblk, name, q, file, line)

    integer, intent(in) :: iblk
    character(*), intent(in) :: name
    real(r8), intent(out), pointer :: q(:,:,:)
    character(*), intent(in) :: file
    integer, intent(in) :: line

    integer idx

    idx = tracer_get_idx(name)
    call tracer_get_array(iblk, idx, q, file, line)

  end subroutine tracer_get_array_name

  subroutine tracer_calc_qm(block)

    type(block_type), intent(in) :: block

    integer i, j, k

    if (.not. allocated(tracers)) return

    associate (mesh => block%mesh, &
               qm   => tracers(block%id)%qm)   ! out
    if (idx_qv > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = tracers(block%id)%q(i,j,k,idx_qv)
          end do
        end do
      end do
    end if
    if (idx_qc > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qc)
          end do
        end do
      end do
    end if
    if (idx_qi > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qi)
          end do
        end do
      end do
    end if
    if (idx_qr > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qr)
          end do
        end do
      end do
    end if
    if (idx_qs > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qs)
          end do
        end do
      end do
    end if
    if (idx_qg > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qg)
          end do
        end do
      end do
    end if
    if (idx_qh > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            qm(i,j,k) = qm(i,j,k) + tracers(block%id)%q(i,j,k,idx_qh)
          end do
        end do
      end do
    end if
    call fill_halo(block%halo, qm, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
    end associate

  end subroutine tracer_calc_qm

  subroutine tracer_get_array_qm(iblk, qm)

    integer, intent(in) :: iblk
    real(r8), intent(out), pointer :: qm(:,:,:)

    associate (mesh => tracers(iblk)%mesh)
    ! NOTE: qm is on mesh. This is different from q which is on filter_mesh.
    qm(mesh%full_ims:mesh%full_ime, &
       mesh%full_jms:mesh%full_jme, &
       mesh%full_kms:mesh%full_kme) => tracers(iblk)%qm
    end associate

  end subroutine tracer_get_array_qm

end module tracer_mod