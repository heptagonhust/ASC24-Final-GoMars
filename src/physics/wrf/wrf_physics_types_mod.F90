module wrf_physics_types_mod

  use const_mod
  use physics_types_mod
  use wrf_namelist_mod
  use rad_rrtmgp_types_mod
  use lsm_noahmp_types_mod

  implicit none

  private

  public wrf_state_type
  public wrf_tend_type

  type, extends(physics_state_type) :: wrf_state_type
    ! Time step index
    integer :: time_step = 0
    ! U wind component at 10 m (m s-1)
    real(r8), allocatable, dimension(:  ) :: u10
    ! V wind component at 10 m (m s-1)
    real(r8), allocatable, dimension(:  ) :: v10
    ! Potential temperature at 2 m (K)
    real(r8), allocatable, dimension(:  ) :: pt2
    ! Specific humidity at 2 m (kg kg-1)
    real(r8), allocatable, dimension(:  ) :: qv2
    ! Upward heat flux at surface (W m-2)
    real(r8), allocatable, dimension(:  ) :: hfx
    ! Upward moisture flux at surface (kg m-2 s-1)
    real(r8), allocatable, dimension(:  ) :: qfx
    ! Bulk Richardson
    real(r8), allocatable, dimension(:  ) :: br
    ! Scalar exchange coefficients
    real(r8), allocatable, dimension(:,:) :: exch_h
    ! Exchange coefficients
    real(r8), allocatable, dimension(:,:) :: exch_m
    ! ???
    real(r8), allocatable, dimension(:  ) :: delta
    ! Sea surface zonal currents (m s-1)
    real(r8), allocatable, dimension(:  ) :: uos
    ! Sea surface meridional currents (m s-1)
    real(r8), allocatable, dimension(:  ) :: vos
    ! Land mask
    real(r8), allocatable, dimension(:  ) :: land
    type(rad_rrtmgp_state_type) rad_rrtmgp
    type(lsm_noahmp_state_type) lsm_noahmp
  contains
    procedure :: init  => wrf_state_init
    procedure :: clear => wrf_state_clear
    final wrf_state_final
  end type wrf_state_type

  type, extends(physics_tend_type) :: wrf_tend_type
    ! Potential temperature tendency due to radition (K s-1)
    real(r8), allocatable, dimension(:,:) :: dptdt_rad
  contains
    procedure :: init  => wrf_tend_init
    procedure :: clear => wrf_tend_clear
    final wrf_tend_final
  end type wrf_tend_type

contains

  subroutine wrf_state_init(this, mesh)

    class(wrf_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%u10   (mesh%ncol          ))
    allocate(this%v10   (mesh%ncol          ))
    allocate(this%pt2   (mesh%ncol          ))
    allocate(this%qv2   (mesh%ncol          ))
    allocate(this%hfx   (mesh%ncol          ))
    allocate(this%qfx   (mesh%ncol          ))
    allocate(this%br    (mesh%ncol          ))
    allocate(this%exch_h(mesh%ncol,mesh%nlev))
    allocate(this%exch_m(mesh%ncol,mesh%nlev))
    allocate(this%delta (mesh%ncol          ))
    allocate(this%uos   (mesh%ncol          ))
    allocate(this%vos   (mesh%ncol          ))
    allocate(this%land  (mesh%ncol          ))

    call this%rad_rrtmgp%init(mesh)
    call this%lsm_noahmp%init(mesh)

    call this%physics_state_init(mesh)

  end subroutine wrf_state_init

  subroutine wrf_state_clear(this)

    class(wrf_state_type), intent(inout) :: this

    this%time_step = 0

    if (allocated(this%u10   )) deallocate(this%u10   )
    if (allocated(this%v10   )) deallocate(this%v10   )
    if (allocated(this%pt2   )) deallocate(this%pt2   )
    if (allocated(this%qv2   )) deallocate(this%qv2   )
    if (allocated(this%hfx   )) deallocate(this%hfx   )
    if (allocated(this%qfx   )) deallocate(this%qfx   )
    if (allocated(this%br    )) deallocate(this%br    )
    if (allocated(this%exch_h)) deallocate(this%exch_h)
    if (allocated(this%exch_m)) deallocate(this%exch_m)
    if (allocated(this%delta )) deallocate(this%delta )
    if (allocated(this%uos   )) deallocate(this%uos   )
    if (allocated(this%vos   )) deallocate(this%vos   )
    if (allocated(this%land  )) deallocate(this%land  )

    call this%lsm_noahmp%clear()

    call this%physics_state_clear()

  end subroutine wrf_state_clear

  subroutine wrf_state_final(this)

    type(wrf_state_type), intent(inout) :: this

    call this%clear()

  end subroutine wrf_state_final

  subroutine wrf_tend_init(this, mesh)

    class(wrf_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    allocate(this%dptdt_rad(mesh%ncol,mesh%nlev))

    call this%physics_tend_init(mesh)

  end subroutine wrf_tend_init

  subroutine wrf_tend_clear(this)

    class(wrf_tend_type), intent(inout) :: this

    if (allocated(this%dptdt_rad)) deallocate(this%dptdt_rad)

    call this%physics_tend_clear()

  end subroutine wrf_tend_clear

  subroutine wrf_tend_final(this)

    type(wrf_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine wrf_tend_final

end module wrf_physics_types_mod
