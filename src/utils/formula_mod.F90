module formula_mod

  use const_mod

  implicit none

  private

  public specific_humidity
  public mixing_ratio
  public potential_temperature
  public temperature
  public virtual_temperature
  public dry_air_density
  public moist_air_density

contains

  pure elemental real(r8) function specific_humidity(qv) result(res)

    real(r8), intent(in) :: qv  ! Mixing ratio

    res = qv / (1 + qv)

  end function specific_humidity

  pure elemental real(r8) function mixing_ratio(sh) result(res)

    real(r8), intent(in) :: sh  ! Specific humidity

    res = sh / (1 - sh)

  end function mixing_ratio

  pure elemental real(r8) function potential_temperature(t, p, qv) result(res)

    real(r8), intent(in) :: t   ! Temperature
    real(r8), intent(in) :: p   ! Pressure
    real(r8), intent(in) :: qv  ! Mixing ratio

    res = t * (p0 / p)**Rd_o_cpd * (1 + Rv_o_Rd * qv)

  end function potential_temperature

  pure elemental real(r8) function temperature(pt, p, qv) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p
    real(r8), intent(in) :: qv

    res = pt * (p / p0)**Rd_o_cpd / (1 + Rv_o_Rd * qv)

  end function temperature

  pure elemental real(r8) function virtual_temperature(t, sh) result(res)

    real(r8), intent(in) :: t   ! Temperature
    real(r8), intent(in) :: sh  ! Specific humidity

    res = t * (1 + (Rv_o_Rd - 1) * sh)

  end function virtual_temperature

  pure elemental real(r8) function dry_air_density(pt, p) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p

    res = p0 / Rd / pt * (p / p0)**cvd_o_cpd

  end function dry_air_density

  pure elemental real(r8) function moist_air_density(t, p, qv) result(res)

    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    real(r8), intent(in) :: qv  ! Mixing ratio

    res = p / Rd / virtual_temperature(t, specific_humidity(qv))

  end function moist_air_density

end module formula_mod
