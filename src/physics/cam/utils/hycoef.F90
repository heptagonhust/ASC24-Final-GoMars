module hycoef

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plevp

  implicit none
  private
  save

  !-----------------------------------------------------------------------
  !
  ! Purpose: Hybrid level definitions: p = a*p0 + b*ps
  !          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
  !          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
  !
  !-----------------------------------------------------------------------

  public hycoef_init
  public hycoef_final

  real(r8), public, allocatable :: hyai(:)
  real(r8), public, allocatable :: hybi(:)
  real(r8), public, allocatable :: hyam(:)
  real(r8), public, allocatable :: hybm(:)
  real(r8), public, allocatable :: hybd(:)         ! difference  in b (hybi) across layers
  real(r8), public, allocatable :: hypi(:)         ! reference pressures at interfaces
  real(r8), public, allocatable :: hypm(:)         ! reference pressures at midpoints
  real(r8), public, allocatable :: hypd(:)         ! reference pressure layer thickness
  real(r8), public, protected :: ps0 = 1.0e5_r8    ! Base state surface pressure (pascals)
  real(r8), public, protected :: psr = 1.0e5_r8    ! Reference surface pressure (pascals)

  integer, public :: nprlev                        ! number of pure pressure levels at top

contains

  subroutine hycoef_init(hyai_in, hybi_in, ps0_in, psr_in)

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Defines the locations of model interfaces from input data in the
    ! hybrid coordinate scheme.  Actual pressure values of model level
    ! interfaces are determined elsewhere from the fields set here.
    !
    ! Method:
    ! the following fields are set:
    ! hyai     fraction of reference pressure used for interface pressures
    ! hyam     fraction of reference pressure used for midpoint pressures
    ! hybi     fraction of surface pressure used for interface pressures
    ! hybm     fraction of surface pressure used for midpoint pressures
    ! hybd     difference of hybi's
    ! hypi     reference state interface pressures
    ! hypm     reference state midpoint pressures
    ! hypd     reference state layer thicknesses
    ! hypdln   reference state layer thicknesses (log p)
    ! hyalph   distance from interface to level (used in integrals)
    ! prsfac   log pressure extrapolation factor (used to compute psl)
    !
    ! Author: B. Boville
    !
    ! Modified by Li Dong at 2023-05-10
    !
    !-----------------------------------------------------------------------

    real(r8), intent(in) :: hyai_in(:)
    real(r8), intent(in) :: hybi_in(:)
    real(r8), intent(in) :: ps0_in
    real(r8), intent(in) :: psr_in

    integer k

    call hycoef_final()

    allocate(hyai(plevp))
    allocate(hybi(plevp))
    allocate(hyam(plev ))
    allocate(hybm(plev ))
    allocate(hybd(plev ))
    allocate(hypi(plevp))
    allocate(hypm(plev ))
    allocate(hypd(plev ))

    hyai = hyai_in
    hybi = hybi_in
    ps0  = ps0_in
    psr  = psr_in
    do k = 1, plev
      hyam(k) = 0.5_r8 * (hyai(k) + hyai(k+1))
      hybm(k) = 0.5_r8 * (hybi(k) + hybi(k+1))
    end do

    nprlev = 0
    do k = 1, plev
      if (nprlev == 0 .and. hybi(k) /= 0.0_r8) then
        nprlev = k - 1
        exit
      end if
    end do

    ! Set nprlev if no nonzero b's have been found. All interfaces are
    ! pure pressure. A pure pressure model requires other changes as well.
    if (nprlev == 0) nprlev = plev + 2

    ! Set delta sigma part of layer thickness and reference state midpoint
    ! pressures
    do k=1,plev
      hybd(k) = hybi(k+1) - hybi(k)
      hypm(k) = hyam(k) * ps0 + hybm(k) * psr
    end do

    ! Reference state interface pressures
    do k = 1, plevp
      hypi(k) = hyai(k) * ps0 + hybi(k) * psr
    end do

    ! Reference state layer thicknesses
    do k = 1, plev
      hypd(k) = hypi(k+1) - hypi(k)
    end do

  end subroutine hycoef_init

  subroutine hycoef_final()

    if (allocated(hyai)) deallocate(hyai)
    if (allocated(hybi)) deallocate(hybi)
    if (allocated(hyam)) deallocate(hyam)
    if (allocated(hybm)) deallocate(hybm)
    if (allocated(hybd)) deallocate(hybd)
    if (allocated(hypi)) deallocate(hypi)
    if (allocated(hypm)) deallocate(hypm)
    if (allocated(hypd)) deallocate(hypd)

  end subroutine hycoef_final

end module hycoef
