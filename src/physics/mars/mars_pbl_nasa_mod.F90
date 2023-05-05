module mars_pbl_nasa_mod

  use const_mod
  use mars_nasa_mod

  implicit none

  private

  public mars_pbl_nasa_init
  public mars_pbl_nasa_run
  public mars_pbl_nasa_final

contains

  subroutine mars_pbl_nasa_init()



  end subroutine mars_pbl_nasa_init

  subroutine mars_pbl_nasa_run(ncol, nlev, rho)

    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    real(r8), intent(in) :: rho(ncol,nlev)

    integer icol, ilev



  end subroutine mars_pbl_nasa_run

  subroutine mars_pbl_nasa_final()

  end subroutine mars_pbl_nasa_final

end module mars_pbl_nasa_mod