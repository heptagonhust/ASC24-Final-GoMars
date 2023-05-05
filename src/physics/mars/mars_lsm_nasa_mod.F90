module mars_lsm_nasa_mod

  use const_mod
  use mars_data_mod

  implicit none

  private

  public mars_lsm_init
  public mars_lsm_run

contains

  subroutine mars_lsm_init()

    integer iblk, icol, ilev

    ! Initialize soil model.
    do iblk = 1, size(mars_data)
      associate (lat      => mars_data(iblk)%lat     , &
                 tin      => mars_data(iblk)%tin     , &
                 soil_rho => mars_data(iblk)%soil_rho, &
                 soil_cp  => mars_data(iblk)%soil_cp , &
                 soil_cn  => mars_data(iblk)%soil_cn )
      ! Set soil density. It would be better to read this from a file.
      do ilev = 1, nlev_soil
        do icol = 1, mars_data(iblk)%ncol
          soil_rho(icol,ilev) = 1481.39_r8
          soil_cp (icol,ilev) = 840.0_r8
        end do
      end do
      do ilev = 1, nlev_soil
        do icol = 1, mars_data(iblk)%ncol
          soil_cn(icol,ilev) = tin(icol,ilev)**2 / (soil_rho(icol,ilev) * soil_cp(icol,ilev))
        end do
      end do
      end associate
    end do

  end subroutine mars_lsm_init

  subroutine mars_lsm_run()

  end subroutine mars_lsm_run

end module mars_lsm_nasa_mod