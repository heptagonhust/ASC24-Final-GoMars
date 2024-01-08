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
  use mars_nasa_spectra_mod

  implicit none

  private

  public mars_nasa_optics_init
  public mars_nasa_optics_final


  ! Absorption coefficient (m-1): βa
  ! Scattering coefficient (m-1): βs
  ! Volume extinction coefficient (m-1): βe = βa + βs
  ! Single scattering albedo (1): ω = βs / βe
  ! Optical depth (1): τ(s1,s2) = ∫ βe(s) ds from s1 to s2
  ! Transmittance (1): t(s1,s2) = exp(-τ(s1,s2))
  ! Mass extinction coefficient: ke = βe / ρ
  ! Scattering asymmetry parameter: g = ∫ p(θ) cos(θ) dΩ / ∫ p(θ) dΩ

  integer, public :: nbin_dst = 0
  integer, public :: nratio_cld = 0
  real(r8), public, allocatable, dimension(    :) :: dst_r_bnds
  real(r8), public, allocatable, dimension(    :) :: dst_be0_vis, dst_be0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_bs0_vis, dst_bs0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_gs0_vis, dst_gs0_ir
  real(r8), public, allocatable, dimension(    :) :: dst_ws0_vis, dst_ws0_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_be_vis , dst_be_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_bs_vis , dst_bs_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_gs_vis , dst_gs_ir
  real(r8), public, allocatable, dimension(  :,:) :: dst_ws_vis , dst_ws_ir
  real(r8), public, allocatable, dimension(:    ) :: cld_ratio
  real(r8), public, allocatable, dimension(:,:,:) :: cld_be_vis , cld_be_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_bs_vis , cld_bs_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_gs_vis , cld_gs_ir
  real(r8), public, allocatable, dimension(:,:,:) :: cld_ws_vis , cld_ws_ir

contains

  subroutine mars_nasa_optics_init(nlev_rad)

    integer, intent(in) :: nlev_rad

    integer i, n

    call mars_nasa_optics_final()

    call fiona_open_dataset('dust_optics', dust_optics_file)
    call fiona_get_dim('dust_optics', 'bin', size=nbin_dst)
    call fiona_get_dim('dust_optics', 'spec_vis', size=n)
    if (n /= spec_vis%n) then
      stop 999
    end if
    call fiona_get_dim('dust_optics', 'spec_ir', size=n)
    if (n /= spec_ir%n) then
      stop 999
    end if
    allocate(dst_r_bnds (nbin_dst+1))
    allocate(dst_be0_vis(         spec_vis%n))
    allocate(dst_bs0_vis(         spec_vis%n))
    allocate(dst_gs0_vis(         spec_vis%n))
    allocate(dst_ws0_vis(         spec_vis%n))
    allocate(dst_be_vis (nbin_dst,spec_vis%n))
    allocate(dst_bs_vis (nbin_dst,spec_vis%n))
    allocate(dst_gs_vis (nbin_dst,spec_vis%n))
    allocate(dst_ws_vis (nbin_dst,spec_vis%n))
    allocate(dst_be0_ir (         spec_ir %n))
    allocate(dst_bs0_ir (         spec_ir %n))
    allocate(dst_gs0_ir (         spec_ir %n))
    allocate(dst_ws0_ir (         spec_ir %n))
    allocate(dst_be_ir  (nbin_dst,spec_ir %n))
    allocate(dst_bs_ir  (nbin_dst,spec_ir %n))
    allocate(dst_gs_ir  (nbin_dst,spec_ir %n))
    allocate(dst_ws_ir  (nbin_dst,spec_ir %n))
    call fiona_start_input('dust_optics')
    call fiona_input('dust_optics', 'dust_r_bnds' , dst_r_bnds )
    call fiona_input('dust_optics', 'dust_be0_vis', dst_be0_vis)
    call fiona_input('dust_optics', 'dust_bs0_vis', dst_bs0_vis)
    call fiona_input('dust_optics', 'dust_gs0_vis', dst_gs0_vis)
    call fiona_input('dust_optics', 'dust_be_vis' , dst_be_vis )
    call fiona_input('dust_optics', 'dust_bs_vis' , dst_bs_vis )
    call fiona_input('dust_optics', 'dust_gs_vis' , dst_gs_vis )
    call fiona_input('dust_optics', 'dust_be0_ir' , dst_be0_ir )
    call fiona_input('dust_optics', 'dust_bs0_ir' , dst_bs0_ir )
    call fiona_input('dust_optics', 'dust_gs0_ir' , dst_gs0_ir )
    call fiona_input('dust_optics', 'dust_be_ir'  , dst_be_ir  )
    call fiona_input('dust_optics', 'dust_bs_ir'  , dst_bs_ir  )
    call fiona_input('dust_optics', 'dust_gs_ir'  , dst_gs_ir  )
    call fiona_end_input('dust_optics')

    do i = 1, spec_vis%n
      if (dst_bs0_vis(i) >= dst_be0_vis(i)) then
        dst_bs0_vis(i) = 0.99999_r8 * dst_be0_vis(i)
      end if
      dst_ws0_vis(i) = dst_bs0_vis(i) / dst_be0_vis(i)
    end do
    do i = 1, spec_ir%n
      if (dst_bs0_ir(i) >= dst_be0_ir(i)) then
        dst_bs0_ir(i) = 0.99999_r8 * dst_be0_ir(i)
      end if
      dst_ws0_ir(i) = dst_bs0_ir(i) / dst_be0_ir(i)
    end do

    call fiona_open_dataset('cld_optics', cld_optics_file)
    call fiona_get_dim('cld_optics', 'ratio', size=nratio_cld)
    call fiona_get_dim('cld_optics', 'spec_vis', size=n)
    if (n /= spec_vis%n) then
      stop 999
    end if
    call fiona_get_dim('cld_optics', 'spec_ir', size=n)
    if (n /= spec_ir%n) then
      stop 999
    end if
    allocate(cld_ratio (nratio_cld))
    allocate(cld_be_vis(nratio_cld,nbin_dst,spec_vis%n))
    allocate(cld_bs_vis(nratio_cld,nbin_dst,spec_vis%n))
    allocate(cld_gs_vis(nratio_cld,nbin_dst,spec_vis%n))
    allocate(cld_ws_vis(nratio_cld,nbin_dst,spec_vis%n))
    allocate(cld_be_ir (nratio_cld,nbin_dst,spec_ir %n))
    allocate(cld_bs_ir (nratio_cld,nbin_dst,spec_ir %n))
    allocate(cld_gs_ir (nratio_cld,nbin_dst,spec_ir %n))
    allocate(cld_ws_ir (nratio_cld,nbin_dst,spec_ir %n))
    call fiona_start_input('cld_optics')
    call fiona_input('cld_optics', 'cld_ratio' , cld_ratio )
    call fiona_input('cld_optics', 'cld_be_vis', cld_be_vis)
    call fiona_input('cld_optics', 'cld_bs_vis', cld_bs_vis)
    call fiona_input('cld_optics', 'cld_gs_vis', cld_gs_vis)
    call fiona_input('cld_optics', 'cld_be_ir' , cld_be_ir )
    call fiona_input('cld_optics', 'cld_bs_ir' , cld_bs_ir )
    call fiona_input('cld_optics', 'cld_gs_ir' , cld_gs_ir )
    call fiona_end_input('cld_optics')

  end subroutine mars_nasa_optics_init

  subroutine mars_nasa_optics_final()

    if (allocated(dst_r_bnds )) deallocate(dst_r_bnds )
    if (allocated(dst_be0_vis)) deallocate(dst_be0_vis)
    if (allocated(dst_bs0_vis)) deallocate(dst_bs0_vis)
    if (allocated(dst_gs0_vis)) deallocate(dst_gs0_vis)
    if (allocated(dst_ws0_vis)) deallocate(dst_ws0_vis)
    if (allocated(dst_be0_ir )) deallocate(dst_be0_ir )
    if (allocated(dst_bs0_ir )) deallocate(dst_bs0_ir )
    if (allocated(dst_gs0_ir )) deallocate(dst_gs0_ir )
    if (allocated(dst_ws0_ir )) deallocate(dst_ws0_ir )
    if (allocated(dst_be_vis )) deallocate(dst_be_vis )
    if (allocated(dst_bs_vis )) deallocate(dst_bs_vis )
    if (allocated(dst_gs_vis )) deallocate(dst_gs_vis )
    if (allocated(dst_ws_vis )) deallocate(dst_ws_vis )
    if (allocated(dst_be_ir  )) deallocate(dst_be_ir  )
    if (allocated(dst_bs_ir  )) deallocate(dst_bs_ir  )
    if (allocated(dst_gs_ir  )) deallocate(dst_gs_ir  )
    if (allocated(dst_ws_ir  )) deallocate(dst_ws_ir  )
    if (allocated(cld_ratio  )) deallocate(cld_ratio  )
    if (allocated(cld_be_vis )) deallocate(cld_be_vis )
    if (allocated(cld_bs_vis )) deallocate(cld_bs_vis )
    if (allocated(cld_gs_vis )) deallocate(cld_gs_vis )
    if (allocated(cld_ws_vis )) deallocate(cld_ws_vis )
    if (allocated(cld_be_ir  )) deallocate(cld_be_ir  )
    if (allocated(cld_bs_ir  )) deallocate(cld_bs_ir  )
    if (allocated(cld_gs_ir  )) deallocate(cld_gs_ir  )
    if (allocated(cld_ws_ir  )) deallocate(cld_ws_ir  )

  end subroutine mars_nasa_optics_final

end module mars_nasa_optics_mod
