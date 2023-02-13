module formula_mod

  use const_mod

  implicit none

  private

  public specific_humidity
  public mixing_ratio
  public potential_temperature
  public temperature
  public virtual_temperature
  public virtual_temperature_from_modified_potential_temperature
  public virtual_potential_temperature
  public dry_air_density
  public moist_air_density
  public buoyancy_frequency
  public local_richardson_number

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

    res = t * (p0 / p)**rd_o_cpd * (1 + rv_o_rd * qv)

  end function potential_temperature

  pure elemental real(r8) function temperature(pt, p, qv) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p
    real(r8), intent(in) :: qv

    res = pt * (p / p0)**rd_o_cpd / (1 + rv_o_rd * qv)

  end function temperature

  pure elemental real(r8) function virtual_temperature(t, qv) result(res)

    real(r8), intent(in) :: t   ! Temperature
    real(r8), intent(in) :: qv  ! Mixing ratio

    res = t * (1 + rv_o_rd * qv) / (1 + qv)

  end function virtual_temperature

  pure elemental real(r8) function virtual_temperature_from_modified_potential_temperature(pt, pk, qv) result(res)

    real(r8), intent(in) :: pt  ! Modified potential temperature
    real(r8), intent(in) :: pk  ! p**(rd/cpd)
    real(r8), intent(in) :: qv  ! Mixing ratio

    res = pt * pk / pk0 / (1 + qv)

  end function virtual_temperature_from_modified_potential_temperature

  pure elemental real(r8) function virtual_potential_temperature(tv, p) result(res)

    real(r8), intent(in) :: tv  ! Virtual temperature
    real(r8), intent(in) :: p   ! Pressure

    res = tv * (p0 / p)**rd_o_cpd

  end function virtual_potential_temperature

  pure elemental real(r8) function dry_air_density(pt, p) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p

    res = p0 / rd / pt * (p / p0)**cvd_o_cpd

  end function dry_air_density

  pure elemental real(r8) function moist_air_density(t, p, qv) result(res)

    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    real(r8), intent(in) :: qv  ! Mixing ratio

    res = p / rd / virtual_temperature(t, specific_humidity(qv))

  end function moist_air_density

  pure elemental real(r8) function buoyancy_frequency(pt1, pt2, z1, z2) result(res)

    real(r8), intent(in) :: pt1
    real(r8), intent(in) :: pt2
    real(r8), intent(in) :: z1
    real(r8), intent(in) :: z2

    res = g * (pt1 - pt2) / (z1 - z2) / (pt1 + pt2) * 2

  end function buoyancy_frequency

  pure elemental real(r8) function local_richardson_number(N2, z1, z2, u1, u2, v1, v2) result(res)

    real(r8), intent(in) :: N2
    real(r8), intent(in) :: z1
    real(r8), intent(in) :: z2
    real(r8), intent(in) :: u1
    real(r8), intent(in) :: u2
    real(r8), intent(in) :: v1
    real(r8), intent(in) :: v2

    real(r8) s2

    s2 = ((u1 - u2)**2 + (v1 - v2)**2) / (z1 - z2)**2

    res = n2 / (s2 + 1.0e-4_r8)

  end function local_richardson_number

end module formula_mod
