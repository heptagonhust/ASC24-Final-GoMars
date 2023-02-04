module moist_mod

  use namelist_mod
  use block_mod
  use adv_mod
  use parallel_mod

  implicit none

  private

  public moist_init
  public moist_link_state
  public calc_qm

contains

  subroutine moist_init()

    call adv_add_tracer('moist', dt_adv, 'qv', 'water vapor' )
    ! call adv_add_tracer('moist', dt_adv, 'ql', 'cloud liquid')
    ! call adv_add_tracer('moist', dt_adv, 'qi', 'cloud ice'   )
    ! call adv_add_tracer('moist', dt_adv, 'qr', 'rain'        )
    ! call adv_add_tracer('moist', dt_adv, 'qs', 'snow'        )
    ! call adv_add_tracer('moist', dt_adv, 'qg', 'graupels'    )
    ! call adv_add_tracer('moist', dt_adv, 'qh', 'hail'        )

  end subroutine moist_init

  subroutine moist_link_state(block)

    type(block_type), intent(inout), target :: block

    integer itime

    associate (mesh      => block%filter_mesh   , &
               adv_batch => block%adv_batches(1))
    do itime = 1, size(block%dstate)
      block%dstate(itime)%qv(        &
        mesh%full_ims:mesh%full_ime, &
        mesh%full_jms:mesh%full_jme, &
        mesh%full_kms:mesh%full_kme) => adv_batch%q(:,:,:,1,adv_batch%old)
      call calc_qm(block, itime)
    end do
    end associate

  end subroutine moist_link_state

  subroutine calc_qm(block, itime)

    type(block_type), intent(inout), target :: block
    integer, intent(in) :: itime

    real(r8), pointer, dimension(:,:,:) :: qv
    integer i, j, k

    associate (mesh => block%filter_mesh     , &
               qv   => block%dstate(itime)%qv, & ! in
               qm   => block%dstate(itime)%qm)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          qm(i,j,k) = qv(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, qm, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
    end associate

  end subroutine calc_qm

end module moist_mod
