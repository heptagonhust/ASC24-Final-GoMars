module mars_nasa_sfc_mod

  use mars_nasa_const_mod
  use mars_nasa_namelist_mod
  use mars_nasa_physics_types_mod
  use mars_nasa_tracers_mod

  implicit none

contains

  subroutine mars_nasa_sfc_run(static, state)

    type(mars_nasa_static_type), intent(in   ) :: static
    type(mars_nasa_state_type ), intent(inout) :: state

    integer icol
    real(r8) alb

    associate (mesh   => state%mesh  , &
               co2ice => state%co2ice, &
               q_sfc  => state%q_sfc , &
               alb    => state%alb   )
    do icol = 1, state%mesh%ncol
      if (co2ice(icol) > 0) then
        alb(icol) = merge(alb_ice_np, alb_ice_sp, mesh%lat(icol) > 0)
      else if (albedo_feedback .and. q_sfc(icol,idx_m_vap) > ice_thresh_kgm2) then
        alb(icol) = ice_albedo
      else
        alb(icol) = static%alb(icol)
      end if
    end do
    end associate

  end subroutine mars_nasa_sfc_run

end module mars_nasa_sfc_mod