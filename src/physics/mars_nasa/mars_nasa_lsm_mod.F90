! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_nasa_lsm_mod

  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_tracers_mod
  use mars_nasa_objects_mod

  implicit none

  private

  public mars_nasa_lsm_init
  public mars_nasa_lsm_final
  public mars_nasa_lsm_run

  real(r8), allocatable, dimension(:) :: soil_z
  real(r8), allocatable, dimension(:) :: soil_z_lev
  real(r8), allocatable, dimension(:) :: soil_dz
  real(r8), allocatable, dimension(:) :: soil_dz_lev

  real(r8), parameter :: factl = 0.25_r8
  real(r8), parameter :: factm = 1.2_r8
  real(r8), parameter :: skind = 0.06_r8
  ! Ground ice depth in the northern hemisphere (m).
  real(r8), parameter :: gidn  = 0.0545_r8
  ! Ground ice depth in the southern hemisphere (m).
  real(r8), parameter :: gids  = 0.0805_r8

contains

  subroutine mars_nasa_lsm_init()

    integer iblk, icol, k

    call mars_nasa_lsm_final()

    allocate(soil_z     (nlev_soil  ))
    allocate(soil_z_lev (nlev_soil+1))
    allocate(soil_dz    (nlev_soil  ))
    allocate(soil_dz_lev(nlev_soil+1))

    ! Setup soil levels.
    soil_dz(1) = factl * skind
    do k = 2, nlev_soil
      soil_dz(k) = factm * soil_dz(k-1)
    end do
    soil_z_lev(1) = 0
    do k = 2, nlev_soil+1
      soil_z_lev(k) = soil_z_lev(k-1) + soil_dz(k-1)
    end do
    do k = 1, nlev_soil
      soil_z(k) = 0.5_r8 * (soil_z_lev(k) + soil_z_lev(k+1))
    end do
    do k = 2, nlev_soil
      soil_dz_lev(k) = 0.5_r8 * (soil_dz(k-1) + soil_dz(k))
    end do

    ! Initialize soil conductivity.
    do iblk = 1, size(objects)
      associate (mesh          => objects(iblk)%mesh               , &
                 gnd_ice       => objects(iblk)%static%gnd_ice     , &
                 tin           => objects(iblk)%static%tin         , &
                 soil_tin      => objects(iblk)%static%soil_tin    , &
                 soil_rho      => objects(iblk)%state%soil_rho     , &
                 soil_cp       => objects(iblk)%state%soil_cp      , &
                 soil_cond     => objects(iblk)%state%soil_cond    , &
                 soil_cond_lev => objects(iblk)%state%soil_cond_lev)
      do icol = 1, objects(iblk)%mesh%ncol
        ! Set soil parameters.
        if (gnd_ice(icol) > 0.5_r8) then
          if (mesh%lat(icol) > 0) then
            do k = 1, nlev_soil
              soil_tin(icol,k) = 2236.995_r8
            end do
          else
            do k = 1, nlev_soil
              soil_tin(icol,k) = 1100.0_r8
            end do
          end if
          do k = 1, nlev_soil
            soil_rho (icol,k) = 1781.99_r8
            soil_cp  (icol,k) = 1404.09_r8
          end do
        else
          do k = 1, nlev_soil
            soil_rho (icol,k) = 1481.39_r8
            soil_cp  (icol,k) = 840.0_r8
          end do
        end if
        do k = 1, nlev_soil
          soil_cond(icol,k) = tin(icol)**2 / (soil_rho(icol,k) * soil_cp(icol,k))
        end do
        soil_cond_lev(icol,1) = soil_cond(icol,1)
        do k = 2, nlev_soil - 1
          soil_cond_lev(icol,k) = 0.5_r8 * (soil_cond(icol,k-1) + soil_cond(icol,k))
        end do
      end do
      end associate
    end do

  end subroutine mars_nasa_lsm_init

  subroutine mars_nasa_lsm_final()

    if (allocated(soil_z     )) deallocate(soil_z     )
    if (allocated(soil_z_lev )) deallocate(soil_z_lev )
    if (allocated(soil_dz    )) deallocate(soil_dz    )
    if (allocated(soil_dz_lev)) deallocate(soil_dz_lev)

  end subroutine mars_nasa_lsm_final

  subroutine mars_nasa_lsm_run(static, state, tend)

    type(mars_nasa_static_type), intent(in   ) :: static
    type(mars_nasa_state_type ), intent(inout) :: state
    type(mars_nasa_tend_type  ), intent(inout) :: tend

    integer icol
    real(r8) tsat

    associate (mesh     => state%mesh    , &
               ps       => state%ps      , & ! in
               co2ice   => state%co2ice  , & ! inout
               q_sfc    => state%q_sfc   , & ! in
               fdnl     => state%fdnl    , & ! in
               rhouch   => state%rhouch  , & ! in
               heat_sfc => state%heat_sfc, & ! in
               alb      => state%alb     , & ! out
               dpsdt    => tend%dpsdt    )   ! out
    do icol = 1, state%mesh%ncol
      ! Calculate surface albedo.
      if (co2ice(icol) > 0) then
        alb(icol) = merge(alb_ice_np, alb_ice_sp, mesh%lat(icol) > 0)
      else if (albedo_feedback .and. q_sfc(icol,idx_m_vap) > ice_thresh_kgm2) then
        alb(icol) = ice_albedo
      else
        alb(icol) = static%alb(icol)
      end if
      ! Calculate CO2 frost point at this surface pressure (hPa).
      tsat = 3182.48_r8 / (23.3494_r8 - log(ps(icol) * 0.01))
      if (co2ice(icol) <= 0) then ! No CO2 on the ground
        co2ice(icol) = 0
        ! Calculate surface temperature (K) from surface energy balance.
        !
        ! f(Ts) =  Fl↓ + (1 - α)Fs↓ - Fconv - Fcond - εσTs⁴
        !
        ! Fconv = ρ u* cp ch (Ts - T)
        ! Fcond =
        !
        ! f'(T) =  - 4εσT³
        !
      end if
    end do
    end associate
    stop 999

  end subroutine mars_nasa_lsm_run

  subroutine calc_t_sfc(state)

    type(mars_nasa_state_type), intent(inout) :: state



  end subroutine calc_t_sfc

end module mars_nasa_lsm_mod