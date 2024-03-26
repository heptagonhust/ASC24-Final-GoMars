! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module damp_mod

  use const_mod
  use namelist_mod
  use block_mod
  use div_damp_mod
  use vor_damp_mod
  use smag_damp_mod
  use laplace_damp_mod

  implicit none

  private

  public damp_init
  public damp_final
  public damp_run

contains

  subroutine damp_init()

    call div_damp_init()
    call vor_damp_init()
    call smag_damp_init()
    call laplace_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call div_damp_final()
    call vor_damp_final()
    call smag_damp_final()
    call laplace_damp_final()

  end subroutine damp_final

  subroutine damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer j

    integer i, k
    real(r8) c
    associate (mesh => block%mesh, u_lon => dstate%u_lon, v_lat => dstate%v_lat)
    ! This nudging of polar v helps to keep the flow neat around the poles.
    ! NOTE: DO NOT REMOVE IT!
    c = 0.8_r8
    do j = mesh%half_jms, mesh%half_jme
      if (mesh%is_south_pole(j)) then
        u_lon%d(:,j+1,:) = c * u_lon%d(:,j+1,:) + (1 - c) * u_lon%d(:,j+2,:)
        v_lat%d(:,j  ,:) = c * v_lat%d(:,j  ,:) + (1 - c) * v_lat%d(:,j+1,:)
      else if (mesh%is_north_pole(j+1)) then
        u_lon%d(:,j  ,:) = c * u_lon%d(:,j  ,:) + (1 - c) * u_lon%d(:,j-1,:)
        v_lat%d(:,j  ,:) = c * v_lat%d(:,j  ,:) + (1 - c) * v_lat%d(:,j-1,:)
      end if
    end do
    end associate

    if (use_vor_damp) then
      do j = 1, vor_damp_cycles
        call vor_damp_run(block, dstate, dt / vor_damp_cycles)
      end do
    end if
    if (use_div_damp) then
      do j = 1, div_damp_cycles
        call div_damp_run(block, dstate, dt / div_damp_cycles)
      end do
    end if
    if (use_smag_damp) then
      do j = 1, smag_damp_cycles
        call smag_damp_run(block, dstate, dt / smag_damp_cycles)
      end do
    end if

  end subroutine damp_run

end module damp_mod
