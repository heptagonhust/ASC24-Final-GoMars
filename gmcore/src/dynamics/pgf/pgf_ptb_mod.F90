! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module implements horizontal pressure gradient force (PGF) by using
!   reference profiles to get perturbed quantities.
!
! History:
!
!   20230821: Initial creation.
!
! Authours:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module pgf_ptb_mod

  use const_mod
  use block_mod
  use formula_mod
  use vert_coord_mod
  use tracer_mod
  use latlon_field_types_mod
  use latlon_parallel_mod

  implicit none

  private

  public pgf_ptb_init_after_ic
  public pgf_ptb_final
  public pgf_ptb_prepare
  public pgf_ptb_run

  type ref_profile_type
    type(latlon_field3d_type) mg
    type(latlon_field3d_type) gz
    type(latlon_field3d_type) gz_lev
    type(latlon_field3d_type) ad ! 1 / rhod
    type(latlon_field3d_type) mpt
    type(latlon_field3d_type) dmgdx
    type(latlon_field3d_type) dmgdy
  contains
    procedure :: init => ref_profile_init
    procedure :: clear => ref_profile_clear
    final :: ref_profile_final
  end type ref_profile_type

  type(ref_profile_type), allocatable :: ref_profiles(:)

contains

  subroutine pgf_ptb_init_after_ic()

    real(r8), parameter :: a   = 0.09923_r8
    real(r8), parameter :: b   = 247.7874_r8
    real(r8), parameter :: c   = -1.0385_r8
    real(r8), parameter :: ps0 = 1013.0_r8 ! hPa
    real(r8) p, t
    integer iblk, i, j, k

    call pgf_ptb_final()

    allocate(ref_profiles(size(blocks)))
    do iblk = 1, size(blocks)
      associate (mesh  => blocks(iblk)%mesh         , &
                 block => blocks(iblk)              , &
                 mgs   => blocks(iblk)%dstate(1)%mgs, &
                 pro   => ref_profiles(iblk)        )
      call pro%init(mesh, blocks(iblk)%halo)
      do k = mesh%half_kds, mesh%half_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            p = vert_coord_calc_mg_lev(k, mgs%d(i,j)) / 100
            pro%gz_lev%d(i,j,k) = -rd * (a * (p - ps0) + b / (1 + c) * (p**(1 + c) - ps0**(1 + c)))
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            p = vert_coord_calc_mg(k, mgs%d(i,j)) / 100
            t = p * (a + b * p**c)
            pro%mg %d(i,j,k) = p * 100
            pro%gz %d(i,j,k) = -rd * (a * (p - ps0) + b / (1 + c) * (p**(1 + c) - ps0**(1 + c)))
            pro%ad %d(i,j,k) = rd * t / (p * 100)
            pro%mpt%d(i,j,k) = dry_potential_temperature(t, pro%mg%d(i,j,k)) * ( &
              vert_coord_calc_mg_lev(k+1, mgs%d(i,j)) - &
              vert_coord_calc_mg_lev(k  , mgs%d(i,j))   &
            )
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            pro%dmgdx%d(i,j,k) = (pro%mg%d(i+1,j,k) - pro%mg%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            pro%dmgdy%d(i,j,k) = (pro%mg%d(i,j+1,k) - pro%mg%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      end associate
    end do

  end subroutine pgf_ptb_init_after_ic

  subroutine pgf_ptb_final()

    if (allocated(ref_profiles)) deallocate(ref_profiles)

  end subroutine pgf_ptb_final

  subroutine pgf_ptb_prepare(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    associate (mesh   => block%mesh            , &
               p      => dstate%p              , & ! in
               gz     => dstate%gz             , & ! in
               rhod   => dstate%rhod           , & ! in
               p_lev  => dstate%p_lev          , & ! in
               dmg    => dstate%dmg            , & ! in
               pro    => ref_profiles(block%id), & ! in
               p_ptb  => block%aux%p_ptb       , & ! out
               gz_ptb => block%aux%gz_ptb      , & ! out
               dp_ptb => block%aux%dp_ptb      , & ! out
               ad_ptb => block%aux%ad_ptb      )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          p_ptb %d(i,j,k) = p %d(i,j,k) - pro%mg%d(i,j,k)
          gz_ptb%d(i,j,k) = gz%d(i,j,k) - pro%gz%d(i,j,k)
          dp_ptb%d(i,j,k) = (p_lev%d(i,j,k+1) - p_lev%d(i,j,k)) - dmg%d(i,j,k)
          ad_ptb%d(i,j,k) = 1.0_r8 / rhod%d(i,j,k) - pro%ad%d(i,j,k)
        end do
      end do
    end do
    end associate

  end subroutine pgf_ptb_prepare

  subroutine pgf_ptb_run(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    real(r8) L, tmp1, tmp2, tmp3, tmp4, tmp
    integer i, j, k

    associate (mesh   => block%mesh                  , &
               qm     => tracers(block%id)%qm        , & ! in
               dmgdx  => ref_profiles(block%id)%dmgdx, & ! in
               dmgdy  => ref_profiles(block%id)%dmgdy, & ! in
               ad_ptb => block%aux%ad_ptb            , & ! in
               p_ptb  => block%aux%p_ptb             , & ! in
               gz_ptb => block%aux%gz_ptb            , & ! in
               dp_ptb => block%aux%dp_ptb            , & ! in
               dmg    => dstate%dmg                  , & ! in
               rhod   => dstate%rhod                 , & ! in
               gz     => dstate%gz                   , & ! in
               du     => dtend%du                    , & ! out
               dv     => dtend%dv                    )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i+1,j,k))
          tmp1 = 0.5_r8 * (ad_ptb%d(i,j,k) + ad_ptb%d(i+1,j,k)) * dmgdx%d(i,j,k)
          tmp2 = 0.5_r8 * (1.0_r8 / rhod%d(i,j,k) + 1.0_r8 / rhod%d(i+1,j,k)) * &
                 (p_ptb%d(i+1,j,k) - p_ptb%d(i,j,k)) / mesh%de_lon(j)
          tmp3 = (gz_ptb%d(i+1,j,k) - gz_ptb%d(i,j,k)) / mesh%de_lon(j)
          tmp4 = 0.5_r8 * (dp_ptb%d(i,j,k) / dmg%d(i,j,k) + dp_ptb%d(i+1,j,k) / dmg%d(i+1,j,k)) * &
                 (gz%d(i+1,j,k) - gz%d(i,j,k)) / mesh%de_lon(j)
          tmp = -(tmp1 + tmp2 + tmp3 + tmp4) / L
          du%d(i,j,k) = du%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_pgf%d(i,j,k) = tmp
#endif
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          L = 1 + 0.5_r8 * (qm%d(i,j,k) + qm%d(i,j+1,k))
          tmp1 = 0.5_r8 * (ad_ptb%d(i,j,k) + ad_ptb%d(i,j+1,k)) * dmgdy%d(i,j,k)
          tmp2 = 0.5_r8 * (1.0_r8 / rhod%d(i,j,k) + 1.0_r8 / rhod%d(i,j+1,k)) * &
                 (p_ptb%d(i,j+1,k) - p_ptb%d(i,j,k)) / mesh%de_lat(j)
          tmp3 = (gz_ptb%d(i,j+1,k) - gz_ptb%d(i,j,k)) / mesh%de_lat(j)
          tmp4 = 0.5_r8 * (dp_ptb%d(i,j,k) / dmg%d(i,j,k) + dp_ptb%d(i,j+1,k) / dmg%d(i,j+1,k)) * &
                 (gz%d(i,j+1,k) - gz%d(i,j,k)) / mesh%de_lat(j)
          tmp = -(tmp1 + tmp2 + tmp3 + tmp4) / L
          dv%d(i,j,k) = dv%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_pgf%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

  end subroutine pgf_ptb_run

  subroutine ref_profile_init(this, mesh, halo)

    class(ref_profile_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)

    character(field_name_len     ) name
    character(field_long_name_len) long_name
    character(field_units_len    ) units

    call this%clear()

    name      = 'ref_pro_mg'
    long_name = 'Reference dry-air pressure'
    units     = 'Pa'
    call this%mg%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_pro_gz'
    long_name = 'Reference geopotential'
    units     = 'm2 s-2'
    call this%gz%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_pro_gz_lev'
    long_name = 'Reference geopotential on half level'
    units     = 'm2 s-2'
    call this%gz_lev%init(name, long_name, units, 'lev', mesh, halo)

    name      = 'ref_pro_ad'
    long_name = 'Reference inverse dry-air density'
    units     = 'kg-1 m3'
    call this%ad%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_pro_mpt'
    long_name = 'Reference mass weighted potential temperature'
    units     = 'K'
    call this%mpt%init(name, long_name, units, 'cell', mesh, halo)

    name      = 'ref_pro_dmgdx'
    long_name = 'Reference dry-air mass flux zonal gradient'
    units     = 'Pa m-2 s-1'
    call this%dmgdx%init(name, long_name, units, 'lon', mesh, halo)

    name      = 'ref_pro_dmgdy'
    long_name = 'Reference dry-air mass flux meridional gradient'
    units     = 'Pa m-2 s-1'
    call this%dmgdy%init(name, long_name, units, 'lat', mesh, halo)

  end subroutine ref_profile_init

  subroutine ref_profile_clear(this)

    class(ref_profile_type), intent(inout) :: this

    call this%mg    %clear()
    call this%gz    %clear()
    call this%gz_lev%clear()
    call this%ad    %clear()
    call this%mpt   %clear()
    call this%dmgdx %clear()
    call this%dmgdy %clear()

  end subroutine ref_profile_clear

  subroutine ref_profile_final(this)

    type(ref_profile_type), intent(inout) :: this

    call this%clear()

  end subroutine ref_profile_final

end module pgf_ptb_mod
