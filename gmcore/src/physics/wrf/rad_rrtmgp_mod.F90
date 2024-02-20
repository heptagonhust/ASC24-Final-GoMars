module rad_rrtmgp_mod

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
  use mo_rte_config
  use wrf_objects_mod

  implicit none

  private

  public rad_rrtmgp_init
  public rad_rrtmgp_final
  public rad_rrtmgp_run

contains

  subroutine rad_rrtmgp_init()

    integer iblk, ierr

  end subroutine rad_rrtmgp_init

  subroutine rad_rrtmgp_final()

  end subroutine rad_rrtmgp_final

  subroutine rad_rrtmgp_run()

  end subroutine rad_rrtmgp_run

end module rad_rrtmgp_mod
