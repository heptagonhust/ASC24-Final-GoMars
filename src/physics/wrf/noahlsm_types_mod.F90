module noahlsm_types_mod

  use const_mod
  use physics_types_mod
  use wrf_namelist_mod

  implicit none

  private

  public noahlsm_state_type

  type noahlsm_state_type
    ! Roof layer temperature (K)
    real(r8), allocatable, dimension(:,:) :: trb_urb4d
    ! Wall layer temperature (K)
    real(r8), allocatable, dimension(:,:) :: tw1_urb4d
    ! Wall layer temperature (K)
    real(r8), allocatable, dimension(:,:) :: tw2_urb4d
    ! Road layer temperature (K)
    real(r8), allocatable, dimension(:,:) :: tgb_urb4d
    ! Indoor temperature (K)
    real(r8), allocatable, dimension(:,:) :: tlev_urb3d
    ! Indoor specific humidity (kg kg-1)
    real(r8), allocatable, dimension(:,:) :: qlev_urb3d
    ! Window temperature (K)
    real(r8), allocatable, dimension(:,:) :: tw1lev_urb3d
    ! Window temperature (K)
    real(r8), allocatable, dimension(:,:) :: tw2lev_urb3d
    ! Ground temperature below building (K)
    real(r8), allocatable, dimension(:,:) :: tglev_urb3d
    ! Floor temperature (K)
    real(r8), allocatable, dimension(:,:) :: tflev_urb3d
    ! Latent heat flux from air conditioner (W m-2)
    real(r8), allocatable, dimension(:  ) :: lf_ac_urb3d
    ! Sensible heat flux from air conditioner (W m-2)
    real(r8), allocatable, dimension(:  ) :: sf_ac_urb3d
    ! Consumption of air conditioner (W m-2)
    real(r8), allocatable, dimension(:  ) :: cm_ac_urb3d
    ! Latent heat flux from urban ventilation (W m-2)
    real(r8), allocatable, dimension(:  ) :: lfvent_urb3d
    ! Sensible heat flux from urban ventilation (W m-2)
    real(r8), allocatable, dimension(:  ) :: sfvent_urb3d
    ! Sensible heat flux from urban surface window (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfwin1_urb3d
    ! Sensible heat flux from urban surface window (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfwin2_urb3d
    ! Sensible heat flux from urban surface (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfw1_urb3d
    ! Sensible heat flux from urban surface (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfw2_urb3d
    ! Sensible heat flux from urban surface (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfr_urb3d
    ! Sensible heat flux from urban surface (W m-2)
    real(r8), allocatable, dimension(:,:) :: sfg_urb3d
  contains
    procedure :: init  => noahlsm_state_init
    procedure :: clear => noahlsm_state_clear
    final noahlsm_state_final
  end type noahlsm_state_type

contains

  subroutine noahlsm_state_init(this, mesh)

    class(noahlsm_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in) :: mesh

    call this%clear()

    allocate(this%trb_urb4d   (mesh%ncol,urban_map_zrd))
    allocate(this%tw1_urb4d   (mesh%ncol,urban_map_zwd))
    allocate(this%tw2_urb4d   (mesh%ncol,urban_map_zwd))
    allocate(this%tgb_urb4d   (mesh%ncol,urban_map_gd ))
    allocate(this%tlev_urb3d  (mesh%ncol,urban_map_bd ))
    allocate(this%qlev_urb3d  (mesh%ncol,urban_map_bd ))
    allocate(this%tw1lev_urb3d(mesh%ncol,urban_map_wd ))
    allocate(this%tw2lev_urb3d(mesh%ncol,urban_map_wd ))
    allocate(this%tglev_urb3d (mesh%ncol,urban_map_gbd))
    allocate(this%tflev_urb3d (mesh%ncol,urban_map_fbd))
    allocate(this%lf_ac_urb3d (mesh%ncol              ))
    allocate(this%sf_ac_urb3d (mesh%ncol              ))
    allocate(this%cm_ac_urb3d (mesh%ncol              ))
    allocate(this%lfvent_urb3d(mesh%ncol              ))
    allocate(this%sfvent_urb3d(mesh%ncol              ))
    allocate(this%sfwin1_urb3d(mesh%ncol,urban_map_wd ))
    allocate(this%sfwin2_urb3d(mesh%ncol,urban_map_wd ))
    allocate(this%sfw1_urb3d  (mesh%ncol,urban_map_zd ))
    allocate(this%sfw2_urb3d  (mesh%ncol,urban_map_zd ))
    allocate(this%sfr_urb3d   (mesh%ncol,urban_map_zdf))
    allocate(this%sfg_urb3d   (mesh%ncol,num_urban_ndm))

  end subroutine noahlsm_state_init

  subroutine noahlsm_state_clear(this)

    class(noahlsm_state_type), intent(inout) :: this

    if (allocated(this%trb_urb4d   )) deallocate(this%trb_urb4d   )
    if (allocated(this%tw1_urb4d   )) deallocate(this%tw1_urb4d   )
    if (allocated(this%tw2_urb4d   )) deallocate(this%tw2_urb4d   )
    if (allocated(this%tgb_urb4d   )) deallocate(this%tgb_urb4d   )
    if (allocated(this%tlev_urb3d  )) deallocate(this%tlev_urb3d  )
    if (allocated(this%qlev_urb3d  )) deallocate(this%qlev_urb3d  )
    if (allocated(this%tw1lev_urb3d)) deallocate(this%tw1lev_urb3d)
    if (allocated(this%tw2lev_urb3d)) deallocate(this%tw2lev_urb3d)
    if (allocated(this%tglev_urb3d )) deallocate(this%tglev_urb3d )
    if (allocated(this%tflev_urb3d )) deallocate(this%tflev_urb3d )
    if (allocated(this%lf_ac_urb3d )) deallocate(this%lf_ac_urb3d )
    if (allocated(this%sf_ac_urb3d )) deallocate(this%sf_ac_urb3d )
    if (allocated(this%cm_ac_urb3d )) deallocate(this%cm_ac_urb3d )
    if (allocated(this%lfvent_urb3d)) deallocate(this%lfvent_urb3d)
    if (allocated(this%sfvent_urb3d)) deallocate(this%sfvent_urb3d)
    if (allocated(this%sfwin1_urb3d)) deallocate(this%sfwin1_urb3d)
    if (allocated(this%sfwin2_urb3d)) deallocate(this%sfwin2_urb3d)
    if (allocated(this%sfw1_urb3d  )) deallocate(this%sfw1_urb3d  )
    if (allocated(this%sfw2_urb3d  )) deallocate(this%sfw2_urb3d  )
    if (allocated(this%sfr_urb3d   )) deallocate(this%sfr_urb3d   )
    if (allocated(this%sfg_urb3d   )) deallocate(this%sfg_urb3d   )

  end subroutine noahlsm_state_clear

  subroutine noahlsm_state_final(this)

    type(noahlsm_state_type), intent(inout) :: this

    call this%clear()

  end subroutine noahlsm_state_final

end module noahlsm_types_mod
