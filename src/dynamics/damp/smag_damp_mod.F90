module smag_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use block_mod

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
      decay_from_top(k) = exp((k - 1)**2 * log(0.01d0) / k0**2) + 1
    end do
    ! FIXME: Disable the decay for the time being.
    decay_from_top = 1

  end subroutine smag_damp_init

  subroutine smag_damp_final()

    if (allocated(decay_from_top)) deallocate(decay_from_top)

  end subroutine smag_damp_final

  subroutine smag_damp_run(block, dt, dtend, dstate)

    type(block_type), intent(inout) :: block
    real(r8), intent(in) :: dt
    type(dtend_type), intent(inout) :: dtend
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) ls2

    associate (mesh      => block%mesh       , &
               smag_t    => block%aux%smag_t , & ! working array
               smag_s    => block%aux%smag_s , & ! working array
               kmh_lon   => block%aux%kmh_lon, & ! working array
               kmh_lat   => block%aux%kmh_lat, & ! working array
               kmh       => block%aux%kmh    , & ! working array
               dudt      => dtend%du         , & ! working array
               dvdt      => dtend%dv         , & ! working array
               dptdt     => dtend%dpt        , & ! working array
               u         => dstate%u_lon     , & ! inout
               v         => dstate%v_lat     , & ! inout
               pt        => dstate%pt        )   ! inout
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          smag_t(i,j,k) = (                       &
            u(i,j,k) - u(i-1,j,k)                 &
          ) / mesh%de_lon(j) - (                  &
            v(i,j  ,k) * mesh%half_cos_lat(j  ) - &
            v(i,j-1,k) * mesh%half_cos_lat(j-1)   &
          ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(block%halo, smag_t, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          smag_s(i,j,k) = (                       &
            v(i+1,j,k) - v(i,j,k)                 &
          ) / mesh%le_lat(j) + (                  &
            u(i,j+1,k) * mesh%full_cos_lat(j+1) - &
            u(i,j  ,k) * mesh%full_cos_lat(j  )   &
          ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(block%halo, smag_s, full_lon=.false., full_lat=.false., full_lev=.true., east_halo=.false., north_halo=.false.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
        do i = mesh%half_ids, mesh%half_ide
          kmh_lon(i,j,k) = ls2 * sqrt(                           &
            0.5_r8 * (smag_t(i,j,k)**2 + smag_t(i+1,j  ,k)**2) + &
            0.5_r8 * (smag_s(i,j,k)**2 + smag_s(i  ,j-1,k)**2)   &
          )
        end do
      end do
    end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        ls2 = smag_damp_coef / (1 / mesh%le_lat(j)**2 + 1 / mesh%de_lat(j)**2) * decay_from_top(k)
        do i = mesh%full_ids, mesh%full_ide
          kmh_lat(i,j,k) = ls2 * sqrt(                           &
            0.5_r8 * (smag_t(i,j,k)**2 + smag_t(i  ,j+1,k)**2) + &
            0.5_r8 * (smag_s(i,j,k)**2 + smag_s(i-1,j  ,k)**2)   &
          )
        end do
      end do
    end do

    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
    !     ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
    !     do i = mesh%full_ids, mesh%full_ide
    !       kmh(i,j,k) = ls2 * sqrt(                        &
    !         smag_t(i,j,k)**2 + 0.25_r8 * (                &
    !           smag_s(i-1,j-1,k)**2 + smag_s(i-1,j,k)**2 + &
    !           smag_s(i  ,j-1,k)**2 + smag_s(i  ,j,k)**2   &
    !         )                                             &
    !       )
    !     end do
    !   end do
    ! end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          dudt(i,j,k) = kmh_lon(i,j,k) * (                                           &
            (u(i-1,j,k) - 2 * u(i,j,k) + u(i+1,j,k)) / mesh%de_lon(j)**2 +           &
            ((u(i,j+1,k) - u(i,j  ,k)) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) - &
             (u(i,j  ,k) - u(i,j-1,k)) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                                &
          )
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u(i,j,k) = u(i,j,k) + dt * dudt(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        if (j == global_mesh%half_jds .or. j == global_mesh%half_jde) then
          do i = mesh%full_ids, mesh%full_ide
            dvdt(i,j,k) = kmh_lat(i,j,k) * (                                           &
              (v(i-1,j,k) - 2 * v(i,j,k) + v(i+1,j,k)) / mesh%le_lat(j)**2             &
            )
          end do
        else
          do i = mesh%full_ids, mesh%full_ide
            dvdt(i,j,k) = kmh_lat(i,j,k) * (                                           &
              (v(i-1,j,k) - 2 * v(i,j,k) + v(i+1,j,k)) / mesh%le_lat(j)**2 +           &
              ((v(i,j+1,k) - v(i,j  ,k)) / mesh%le_lon(j+1) * mesh%full_cos_lat(j+1) - &
               (v(i,j  ,k) - v(i,j-1,k)) / mesh%le_lon(j  ) * mesh%full_cos_lat(j  )   &
              ) / mesh%de_lat(j) / mesh%half_cos_lat(j)                                &
            )
          end do
        end if
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v(i,j,k) = v(i,j,k) + dt * dvdt(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)

    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
    !     do i = mesh%full_ids, mesh%full_ide
    !       dptdt(i,j,k) = kmh(i,j,k) * (                                                &
    !         (pt(i-1,j,k) - 2 * pt(i,j,k) + pt(i+1,j,k)) / mesh%de_lon(j)**2 +          &
    !         ((pt(i,j+1,k) - pt(i,j  ,k)) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) - &
    !          (pt(i,j  ,k) - pt(i,j-1,k)) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)   &
    !         ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                                  &
    !       )
    !     end do
    !   end do
    ! end do
    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
    !     do i = mesh%full_ids, mesh%full_ide
    !       pt(i,j,k) = pt(i,j,k) + dt * dptdt(i,j,k)
    !     end do
    !   end do
    ! end do
    ! call fill_halo(block%halo, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine smag_damp_run

end module smag_damp_mod
