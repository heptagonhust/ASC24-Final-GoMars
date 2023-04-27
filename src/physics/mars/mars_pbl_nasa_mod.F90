module mars_pbl_nasa_mod

  use const_mod
  use mars_data_nasa_mod

  implicit none

  private

  public mars_pbl_nasa_init
  public mars_pbl_nasa_run

  real(r8), parameter :: ric_pbl    = 0.195_r8    ! Critical Richardson number for PBL

contains

  subroutine mars_pbl_nasa_init()



  end subroutine mars_pbl_nasa_init

  subroutine mars_pbl_nasa_run(ncol, nlev, rho)

    integer, intent(in) :: ncol
    integer, intent(in) :: nlev
    real(r8), intent(in) :: rho(ncol,nlev)

    integer icol, ilev



  end subroutine mars_pbl_nasa_run

end module mars_pbl_nasa_mod