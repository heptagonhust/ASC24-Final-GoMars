module clm_comp_shr

  use esmf
  use shr_kind_mod, only: cl => shr_kind_cl

  implicit none

  type(ESMF_Mesh), pointer :: mesh => null()
  type(ESMF_Clock) model_clock
  character(cl) model_meshfile

end module clm_comp_shr
