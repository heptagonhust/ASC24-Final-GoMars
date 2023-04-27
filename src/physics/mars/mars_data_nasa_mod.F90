module mars_data_nasa_mod

  use const_mod

  implicit none

  type mars_data_type
    real(r8), allocatable, dimension(:,:) :: rho_soil   ! Soil density (kg m-3)
    real(r8), allocatable, dimension(:,:) :: c_soil     ! Soil specific capacity (J kg-1 K-1)
  end type mars_data_type

  type(mars_data_type), allocatable :: mars_data(:)

contains

  subroutine mars_data_nasa_init(nblk, ncol, nlev)

    integer, intent(in) :: nblk
    integer, intent(in) :: ncol(:)
    integer, intent(in) :: nlev

    integer iblk

    call mars_data_nasa_final()

    allocate(mars_data(nblk))

    do iblk = 1, nblk
      allocate(mars_data(iblk)%rho_soil(ncol(iblk),nlev))
      allocate(mars_data(iblk)%c_soil  (ncol(iblk),nlev))
    end do

  end subroutine mars_data_nasa_init

  subroutine mars_data_nasa_final()

    integer iblk

    if (allocated(mars_data)) then
      do iblk = 1, size(mars_data)
        if (allocated(mars_data(iblk)%rho_soil)) deallocate(mars_data(iblk)%rho_soil)
        if (allocated(mars_data(iblk)%c_soil  )) deallocate(mars_data(iblk)%c_soil  )
      end do
      deallocate(mars_data)
    end if

  end subroutine mars_data_nasa_final

end module mars_data_nasa_mod