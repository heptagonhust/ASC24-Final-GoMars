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
  public optics_dst_vis
  public optics_dst_ir
  public optics_cld_vis
  public optics_cld_ir


  ! Absorption coefficient (m-1): βa
  ! Scattering coefficient (m-1): βs
  ! Volume extinction coefficient (m-1): βe = βa + βs
  ! Single scattering albedo (1): ω = βs / βe
  ! Optical depth (1): τ(s1,s2) = ∫ βe(s) ds from s1 to s2
  ! Transmittance (1): t(s1,s2) = exp(-τ(s1,s2))
  ! Mass extinction coefficient: ke = βe / ρ

  type optics_type
    ! Particle radius bin number
    integer :: nbin  = 0
    ! Spectral number
    integer :: nspec = 0
    ! Particle radius bin boundaries
    real(r8), allocatable, dimension(:  ) :: r_bnds
    ! Extinction efficient at each spectral interval
    real(r8), allocatable, dimension(  :) :: be0
    ! Scattering efficient at each spectral interval
    real(r8), allocatable, dimension(  :) :: bs0
    ! Asymmetry parameter at each spectral interval
    real(r8), allocatable, dimension(  :) :: gs0
    ! Single scattering albedo at each spectral interval
    real(r8), allocatable, dimension(  :) :: ws0
    ! Extinction efficient at each spectral interval and each particle bin
    real(r8), allocatable, dimension(:,:) :: be
    ! Scattering efficient at each spectral interval and each particle bin
    real(r8), allocatable, dimension(:,:) :: bs
    ! Asymmetry parameter at each spectral interval and each particle bin
    real(r8), allocatable, dimension(:,:) :: gs
  contains
    procedure :: init  => optics_init
    procedure :: clear => optics_clear
    final optics_final
  end type optics_type

  type(optics_type) optics_dst_vis
  type(optics_type) optics_dst_ir
  type(optics_type) optics_cld_vis
  type(optics_type) optics_cld_ir

contains

  subroutine mars_nasa_optics_init(nlev_rad)

    integer, intent(in) :: nlev_rad

    call optics_dst_vis%init( dust_optics_file, 'dust' , 'vis')
    call optics_dst_ir %init( dust_optics_file, 'dust' , 'ir' )
    ! call optics_cld_vis%init(cloud_optics_file, 'cloud', 'vis')
    ! call optics_cld_ir %init(cloud_optics_file, 'cloud', 'ir' )

  end subroutine mars_nasa_optics_init

  subroutine mars_nasa_optics_final()

    call optics_dst_vis%clear()
    call optics_dst_ir %clear()
    call optics_cld_vis%clear()
    call optics_cld_ir %clear()

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
    allocate(this%r_bnds(this%nbin+1))
    allocate(this%be0(          this%nspec))
    allocate(this%bs0(          this%nspec))
    allocate(this%gs0(          this%nspec))
    allocate(this%ws0(          this%nspec))
    allocate(this%be (this%nbin,this%nspec))
    allocate(this%bs (this%nbin,this%nspec))
    allocate(this%gs (this%nbin,this%nspec))
    call fiona_start_input(tag)
    call fiona_input(tag, trim(species) // '_r_bnds', this%r_bnds)
    call fiona_input(tag, trim(species) // '_be0_'// trim(spectra), this%be0)
    call fiona_input(tag, trim(species) // '_bs0_'// trim(spectra), this%bs0)
    call fiona_input(tag, trim(species) // '_gs0_'// trim(spectra), this%gs0)
    call fiona_input(tag, trim(species) // '_be_' // trim(spectra), this%be )
    call fiona_input(tag, trim(species) // '_bs_' // trim(spectra), this%bs )
    call fiona_input(tag, trim(species) // '_gs_' // trim(spectra), this%gs )
    call fiona_end_input(tag)

    do i = 1, this%nspec
      if (this%bs0(i) >= this%be0(i)) then
        this%bs0(i) = 0.99999_r8 * this%be0(i)
      end if
      this%ws0(i) = this%bs0(i) / this%be0(i)
    end do

  end subroutine optics_init

  subroutine optics_clear(this)

    class(optics_type), intent(inout) :: this

    this%nbin  = 0
    this%nspec = 0
    if (allocated(this%r_bnds)) deallocate(this%r_bnds)
    if (allocated(this%be0   )) deallocate(this%be0   )
    if (allocated(this%bs0   )) deallocate(this%bs0   )
    if (allocated(this%gs0   )) deallocate(this%gs0   )
    if (allocated(this%ws0   )) deallocate(this%ws0   )
    if (allocated(this%be    )) deallocate(this%be    )
    if (allocated(this%bs    )) deallocate(this%bs    )
    if (allocated(this%gs    )) deallocate(this%gs    )

  end subroutine optics_clear

  subroutine optics_final(this)

    type(optics_type), intent(inout) :: this

    call this%clear()

  end subroutine optics_final

end module mars_nasa_optics_mod
