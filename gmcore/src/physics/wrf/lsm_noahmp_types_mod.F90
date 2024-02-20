module lsm_noahmp_types_mod

  use const_mod
  use physics_types_mod
  use wrf_namelist_mod
  use NoahmpIOVarType
  use NoahmpIOVarInitMod

  implicit none

  private

  public lsm_noahmp_state_type

  type, extends(NoahmpIO_type) :: lsm_noahmp_state_type
  contains
    procedure :: init  => lsm_noahmp_state_init
    procedure :: clear => lsm_noahmp_state_clear
    final lsm_noahmp_state_final
  end type lsm_noahmp_state_type

contains

  subroutine lsm_noahmp_state_init(this, mesh)

    class(lsm_noahmp_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    this%xstart         = 1
    this%xend           = mesh%ncol
    this%ystart         = 1
    this%yend           = 1
    this%kds            = 1
    this%kde            = mesh%nlev
    this%nsoil          = num_soil_layers
    this%nsnow          = num_snow_layers
    this%urban_map_fbd  = urban_map_fbd
    this%urban_map_gbd  = urban_map_gbd
    this%urban_map_zdf  = urban_map_zdf
    this%urban_map_zrd  = urban_map_zrd
    this%urban_map_zwd  = urban_map_zwd
    this%urban_map_zgrd = urban_map_zgrd
    this%urban_map_bd   = urban_map_bd
    this%urban_map_gd   = urban_map_gd
    this%urban_map_wd   = urban_map_wd
    this%urban_map_zd   = urban_map_zd
    this%num_urban_ndm  = num_urban_ndm
    this%num_urban_hi   = num_urban_hi

    call NoahmpIOVarInitDefault(this)

  end subroutine lsm_noahmp_state_init

  subroutine lsm_noahmp_state_clear(this)

    class(lsm_noahmp_state_type), intent(inout) :: this

  end subroutine lsm_noahmp_state_clear

  subroutine lsm_noahmp_state_final(this)

    type(lsm_noahmp_state_type), intent(inout) :: this

    call this%clear()

  end subroutine lsm_noahmp_state_final

end module lsm_noahmp_types_mod
