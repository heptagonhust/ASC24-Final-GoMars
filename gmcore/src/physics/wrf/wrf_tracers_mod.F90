module wrf_tracers_mod

  use const_mod
  use tracer_mod

  implicit none

  private

  public wrf_tracers_init

contains

  subroutine wrf_tracers_init(dt)

    real(r8), intent(in) :: dt

    call tracer_add('wrf', dt, 'qv', 'Water vapor', 'kg kg-1')

  end subroutine wrf_tracers_init

end module wrf_tracers_mod
