module rad_rrtmgp_types_mod

  use mo_rte_kind
  use mo_optical_props
  use mo_gas_optics_rrtmgp
  use mo_cloud_optics_rrtmgp
  use mo_aerosol_optics_rrtmgp_merra
  use mo_gas_concentrations
  use mo_source_functions
  use mo_fluxes
  use mo_rte_lw
  use mo_rte_sw
  use physics_mesh_mod

  implicit none

  private

  public rad_rrtmgp_state_type

  type rad_rrtmgp_state_type
    type(ty_source_func_lw) lw_src
    type(ty_gas_optics_rrtmgp) k_dist
    type(ty_cloud_optics_rrtmgp) cld_optics
    type(ty_gas_concs) gas_concs
  contains
    procedure init  => rad_rrtmgp_state_init
    procedure clear => rad_rrtmgp_state_clear
    final rad_rrtmgp_state_final
  end type rad_rrtmgp_state_type

contains

  subroutine rad_rrtmgp_state_init(this, mesh)

    class(rad_rrtmgp_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in) :: mesh

  end subroutine rad_rrtmgp_state_init

  subroutine rad_rrtmgp_state_clear(this)

    class(rad_rrtmgp_state_type), intent(inout) :: this

  end subroutine rad_rrtmgp_state_clear

  subroutine rad_rrtmgp_state_final(this)

    type(rad_rrtmgp_state_type), intent(inout) :: this

  end subroutine rad_rrtmgp_state_final

end module rad_rrtmgp_types_mod
