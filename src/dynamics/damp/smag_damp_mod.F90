! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module smag_damp_mod

  use flogger
  use string
  use const_mod
  use math_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod
  use tracer_mod

  implicit none

  private

  public smag_damp_init
  public smag_damp_run
  public smag_damp_final

  real(r8), allocatable, dimension(:), target :: decay_from_top

contains

  subroutine smag_damp_init()

    integer k, k0

    call smag_damp_final()

    allocate(decay_from_top(global_mesh%full_nlev))

    k0 = 8
    do k = global_mesh%full_kds, global_mesh%full_kde
      decay_from_top(k) = exp_two_values(5.0_r8, 1.0_r8, 1.0_r8, real(k0, r8), real(k, r8))
    end do
    ! FIXME: Disable the decay for the time being.
    decay_from_top = 1

  end subroutine smag_damp_init

  subroutine smag_damp_final()

    if (allocated(decay_from_top)) deallocate(decay_from_top)

  end subroutine smag_damp_final

  subroutine smag_damp_run(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k, m
    real(r8) ls2
    real(r8), pointer :: q(:,:,:,:)

    associate (mesh      => block%mesh          , &
               smag_t    => block%aux%smag_t    , & ! working array
               smag_s    => block%aux%smag_s    , & ! working array
               kmh_lon   => block%aux%kmh_lon   , & ! working array
               kmh_lat   => block%aux%kmh_lat   , & ! working array
               kmh       => block%aux%kmh       , & ! working array
               dmg       => dstate%dmg          , & ! working array
               u         => dstate%u_lon        , & ! inout
               v         => dstate%v_lat        , & ! inout
               pt        => dstate%pt           )   ! inout
    ! Horizontal tension strain on centers
    ! ∂u   ∂v
    ! -- - --
    ! ∂x   ∂y
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          smag_t%d(i,j,k) = (                       &
            u%d(i,j,k) - u%d(i-1,j,k)               &
          ) / mesh%de_lon(j) - (                    &
            v%d(i,j  ,k) * mesh%half_cos_lat(j  ) - &
            v%d(i,j-1,k) * mesh%half_cos_lat(j-1)   &
          ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(smag_t, west_halo=.false., south_halo=.false.)

    ! Horizontal shearing strain on vertices
    ! ∂u   ∂v
    ! -- + --
    ! ∂y   ∂x
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          smag_s%d(i,j,k) = (                       &
            v%d(i+1,j,k) - v%d(i,j,k)               &
          ) / mesh%le_lat(j) + (                    &
            u%d(i,j+1,k) * mesh%full_cos_lat(j+1) - &
            u%d(i,j  ,k) * mesh%full_cos_lat(j  )   &
          ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(smag_s, east_halo=.false., north_halo=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
        do i = mesh%half_ids, mesh%half_ide
          kmh_lon%d(i,j,k) = ls2 * sqrt(                             &
            0.5_r8 * (smag_t%d(i,j,k)**2 + smag_t%d(i+1,j  ,k)**2) + &
            0.5_r8 * (smag_s%d(i,j,k)**2 + smag_s%d(i  ,j-1,k)**2)   &
          )
        end do
      end do
    end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        ls2 = smag_damp_coef / (1 / mesh%le_lat(j)**2 + 1 / mesh%de_lat(j)**2) * decay_from_top(k)
        do i = mesh%full_ids, mesh%full_ide
          kmh_lat%d(i,j,k) = ls2 * sqrt(                             &
            0.5_r8 * (smag_t%d(i,j,k)**2 + smag_t%d(i  ,j+1,k)**2) + &
            0.5_r8 * (smag_s%d(i,j,k)**2 + smag_s%d(i-1,j  ,k)**2)   &
          )
        end do
      end do
    end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
        do i = mesh%full_ids, mesh%full_ide
          kmh%d(i,j,k) = ls2 * sqrt(                          &
            smag_t%d(i,j,k)**2 + 0.25_r8 * (                  &
              smag_s%d(i-1,j-1,k)**2 + smag_s%d(i-1,j,k)**2 + &
              smag_s%d(i  ,j-1,k)**2 + smag_s%d(i  ,j,k)**2   &
            )                                                 &
          )
        end do
      end do
    end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u%d(i,j,k) = u%d(i,j,k) + kmh_lon%d(i,j,k) * (                                 &
            (u%d(i-1,j,k) - 2 * u%d(i,j,k) + u%d(i+1,j,k)) / mesh%de_lon(j)**2 +         &
            ((u%d(i,j+1,k) - u%d(i,j  ,k)) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) - &
             (u%d(i,j  ,k) - u%d(i,j-1,k)) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                                    &
          ) * dt
        end do
      end do
    end do
    call fill_halo(u)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        if (j == global_mesh%half_jds .or. j == global_mesh%half_jde) then
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = v%d(i,j,k) + kmh_lat%d(i,j,k) * (                       &
              (v%d(i-1,j,k) - 2 * v%d(i,j,k) + v%d(i+1,j,k)) / mesh%le_lat(j)**2 &
            ) * dt
          end do
        else
          do i = mesh%full_ids, mesh%full_ide
            v%d(i,j,k) = v%d(i,j,k) + kmh_lat%d(i,j,k) * (                                 &
              (v%d(i-1,j,k) - 2 * v%d(i,j,k) + v%d(i+1,j,k)) / mesh%le_lat(j)**2 +         &
              ((v%d(i,j+1,k) - v%d(i,j  ,k)) / mesh%le_lon(j+1) * mesh%full_cos_lat(j+1) - &
               (v%d(i,j  ,k) - v%d(i,j-1,k)) / mesh%le_lon(j  ) * mesh%full_cos_lat(j  )   &
              ) / mesh%de_lat(j) / mesh%half_cos_lat(j)                                    &
            ) * dt
          end do
        end if
      end do
    end do
    call fill_halo(v)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          pt%d(i,j,k) = pt%d(i,j,k) + kmh%d(i,j,k) * (                    &
            (dmg%d(i-1,j,k) * pt%d(i-1,j,k) - 2 *                         &
             dmg%d(i  ,j,k) * pt%d(i  ,j,k) +                             &
             dmg%d(i+1,j,k) * pt%d(i+1,j,k)) / mesh%de_lon(j)**2 +        &
            ((dmg%d(i,j+1,k) * pt%d(i,j+1,k) - dmg%d(i,j,k) * pt%d(i,j,k) &
             ) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) -              &
             (dmg%d(i,j,k) * pt%d(i,j,k) - dmg%d(i,j-1,k) * pt%d(i,j-1,k) &
             ) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)                &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                     &
          ) / dmg%d(i,j,k) * dt
        end do
      end do
    end do
    call fill_halo(pt)
    end associate

  end subroutine smag_damp_run

end module smag_damp_mod
