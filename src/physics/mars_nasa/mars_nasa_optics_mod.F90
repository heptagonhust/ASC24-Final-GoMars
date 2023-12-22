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

module mars_nasa_optics_mod

  use fiona
  use const_mod
  use mars_nasa_namelist_mod

  implicit none

  private

  public mars_nasa_optics_init
  public mars_nasa_optics_final
  public dust_vis
  public dust_ir

  type optics_type
    integer :: nbin  = 0
    integer :: nspec = 0
    real(r8), allocatable, dimension(  :) :: qext0
    real(r8), allocatable, dimension(  :) :: qscat0
    real(r8), allocatable, dimension(  :) :: gscat0
    real(r8), allocatable, dimension(:,:) :: qext_bin
    real(r8), allocatable, dimension(:,:) :: qscat_bin
    real(r8), allocatable, dimension(:,:) :: gscat_bin
    real(r8), allocatable, dimension(  :) :: qext       ! Extinction efficiency at each spectral interval
    real(r8), allocatable, dimension(  :) :: qscat      ! Scattering efficiency at each spectral interval
    real(r8), allocatable, dimension(  :) :: gscat      ! Asymmetry parameter at each spectral interval
    real(r8), allocatable, dimension(  :) :: wscat      ! Single scattering albedo at each spectral interval
  contains
    procedure :: init  => optics_init
    procedure :: clear => optics_clear
    final optics_final
  end type optics_type

  type(optics_type) dust_vis
  type(optics_type) dust_ir

contains

  subroutine mars_nasa_optics_init()

    call dust_vis%init(dust_optics_file, 'dust', 'vis')
    call dust_ir %init(dust_optics_file, 'dust', 'ir' )

  end subroutine mars_nasa_optics_init

  subroutine mars_nasa_optics_final()

    call dust_vis%clear()
    call dust_ir %clear()

  end subroutine mars_nasa_optics_final

  subroutine optics_init(this, optics_file, species, spectra)

    class(optics_type), intent(inout) :: this
    character(*), intent(in) :: optics_file
    character(*), intent(in) :: species
    character(*), intent(in) :: spectra

    character(30) tag
    integer i

    tag = 'optics_' // trim(species) // '_' // trim(spectra)
    call fiona_open_dataset(tag, file_path=optics_file)
    call fiona_get_dim(tag, 'bin', size=this%nbin)
    call fiona_get_dim(tag, 'spec_' // trim(spectra), size=this%nspec)
    allocate(this%qext0    (          this%nspec))
    allocate(this%qscat0   (          this%nspec))
    allocate(this%gscat0   (          this%nspec))
    allocate(this%qext_bin (this%nbin,this%nspec))
    allocate(this%qscat_bin(this%nbin,this%nspec))
    allocate(this%gscat_bin(this%nbin,this%nspec))
    allocate(this%qext     (          this%nspec))
    allocate(this%qscat    (          this%nspec))
    allocate(this%gscat    (          this%nspec))
    allocate(this%wscat    (          this%nspec))
    call fiona_start_input(tag)
    call fiona_input(tag, trim(species) // '_qext0_'     // trim(spectra), this%qext0    )
    call fiona_input(tag, trim(species) // '_qscat0_'    // trim(spectra), this%qscat0   )
    call fiona_input(tag, trim(species) // '_gscat0_'    // trim(spectra), this%gscat0   )
    call fiona_input(tag, trim(species) // '_qext_bin_'  // trim(spectra), this%qext_bin )
    call fiona_input(tag, trim(species) // '_qscat_bin_' // trim(spectra), this%qscat_bin)
    call fiona_input(tag, trim(species) // '_gscat_bin_' // trim(spectra), this%gscat_bin)
    call fiona_end_input(tag)

    do i = 1, this%nspec
      this%qext (i) = this%qext0 (i)
      this%qscat(i) = this%qscat0(i)
      if (this%qscat(i) >= this%qext(i)) then
        this%qscat(i) = 0.99999_r8 * this%qext(i)
      end if
      this%wscat(i) = this%qscat(i) / this%qext(i)
      this%gscat(i) = this%gscat0(i)
    end do

  end subroutine optics_init

  subroutine optics_clear(this)

    class(optics_type), intent(inout) :: this

    if (allocated(this%qext0    )) deallocate(this%qext0    )
    if (allocated(this%qscat0   )) deallocate(this%qscat0   )
    if (allocated(this%gscat0   )) deallocate(this%gscat0   )
    if (allocated(this%qext_bin )) deallocate(this%qext_bin )
    if (allocated(this%qscat_bin)) deallocate(this%qscat_bin)
    if (allocated(this%gscat_bin)) deallocate(this%gscat_bin)
    if (allocated(this%qext     )) deallocate(this%qext     )
    if (allocated(this%qscat    )) deallocate(this%qscat    )
    if (allocated(this%gscat    )) deallocate(this%gscat    )
    if (allocated(this%wscat    )) deallocate(this%wscat    )

  end subroutine optics_clear

  subroutine optics_final(this)

    type(optics_type), intent(inout) :: this

    call this%clear()

  end subroutine optics_final

end module mars_nasa_optics_mod
