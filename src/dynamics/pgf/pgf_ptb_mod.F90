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
  use allocator_mod
  use vert_coord_mod
  use latlon_parallel_mod

  implicit none

  private

  public pgf_ptb_init_after_ic
  public pgf_ptb_final
  public pgf_ptb_prepare
  public pgf_ptb_run

  type ref_profile_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:,:) :: mg
    real(r8), allocatable, dimension(:,:,:) :: gz
    real(r8), allocatable, dimension(:,:,:) :: gz_lev
    real(r8), allocatable, dimension(:,:,:) :: alp
    real(r8), allocatable, dimension(:,:,:) :: mpt
    real(r8), allocatable, dimension(:,:,:) :: dmgdx
    real(r8), allocatable, dimension(:,:,:) :: dmgdy
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
      call pro%init(mesh)
      do k = mesh%half_kds, mesh%half_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            p = vert_coord_calc_mg_lev(k, mgs(i,j)) / 100
            pro%gz_lev(i,j,k) = -rd * (a * (p - ps0) + b / (1 + c) * (p**(1 + c) - ps0**(1 + c)))
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            p = vert_coord_calc_mg(k, mgs(i,j)) / 100
            t = p * (a + b * p**c)
            pro%mg (i,j,k) = p * 100
            pro%gz (i,j,k) = -rd * (a * (p - ps0) + b / (1 + c) * (p**(1 + c) - ps0**(1 + c)))
            pro%alp(i,j,k) = rd * t / (p * 100)
            pro%mpt(i,j,k) = dry_potential_temperature(t, p) * ( &
              vert_coord_calc_mg_lev(k+1, mgs(i,j)) - &
              vert_coord_calc_mg_lev(k  , mgs(i,j))   &
            )
          end do
        end do
      end do
      call fill_halo(block%halo, pro%mg , full_lon=.true., full_lat=.true., full_lev=.true.)
      call fill_halo(block%halo, pro%alp, full_lon=.true., full_lat=.true., full_lev=.true.)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            pro%dmgdx(i,j,k) = (pro%mg(i+1,j,k) - pro%mg(i,j,k)) / mesh%de_lon(j)
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            pro%dmgdy(i,j,k) = (pro%mg(i,j+1,k) - pro%mg(i,j,k)) / mesh%de_lat(j)
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

    associate (mesh    => block%mesh            , &
               p       => dstate%p              , &
               gz      => dstate%gz             , &
               rhod    => dstate%rhod           , &
               p_lev   => dstate%p_lev          , &
               dmg     => dstate%dmg            , &
               pro     => ref_profiles(block%id), &
               p_ptb   => block%aux%p_ptb       , &
               gz_ptb  => block%aux%gz_ptb      , &
               dp_ptb  => block%aux%dp_ptb      , &
               alp_ptb => block%aux%alp_ptb     )
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          p_ptb  (i,j,k) = p (i,j,k) - pro%mg(i,j,k)
          gz_ptb (i,j,k) = gz(i,j,k) - pro%gz(i,j,k)
          dp_ptb (i,j,k) = (p_lev(i,j,k+1) - p_lev(i,j,k)) - dmg(i,j,k)
          alp_ptb(i,j,k) = 1.0_r8 / rhod(i,j,k) - pro%alp(i,j,k)
        end do
      end do
    end do
    call fill_halo(block%halo, p_ptb  , full_lon=.true., full_lat=.true., full_lev=.true., &
                   west_halo=.false., south_halo=.false.)
    call fill_halo(block%halo, gz_ptb , full_lon=.true., full_lat=.true., full_lev=.true., &
                   west_halo=.false., south_halo=.false.)
    call fill_halo(block%halo, dp_ptb , full_lon=.true., full_lat=.true., full_lev=.true., &
                   west_halo=.false., south_halo=.false.)
    call fill_halo(block%halo, alp_ptb, full_lon=.true., full_lat=.true., full_lev=.true., &
                   west_halo=.false., south_halo=.false.)
    end associate

  end subroutine pgf_ptb_prepare

  subroutine pgf_ptb_run(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    real(r8) tmp1, tmp2, tmp3, tmp4
    integer i, j, k

    associate (mesh    => block%mesh                  , &
               alp_ref => ref_profiles(block%id)%alp  , & ! in
               dmgdx   => ref_profiles(block%id)%dmgdx, & ! in
               dmgdy   => ref_profiles(block%id)%dmgdy, & ! in
               alp_ptb => block%aux%alp_ptb           , & ! in
               p_ptb   => block%aux%p_ptb             , & ! in
               gz_ptb  => block%aux%gz_ptb            , & ! in
               dp_ptb  => block%aux%dp_ptb            , & ! in
               dmg     => dstate%dmg                  , & ! in
               gz      => dstate%gz                   , & ! in
               du      => dtend%du                    , & ! out
               dv      => dtend%dv                    )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp1 = 0.5_r8 * (alp_ptb(i,j,k) + alp_ptb(i+1,j,k)) * dmgdx(i,j,k)
          tmp2 = 0.5_r8 * (alp_ref(i,j,k) + alp_ref(i+1,j,k)) * &
                 (p_ptb(i+1,j,k) - p_ptb(i,j,k)) / mesh%de_lon(j)
          tmp3 = (gz_ptb(i+1,j,k) - gz_ptb(i,j,k)) / mesh%de_lon(j)
          tmp4 = 0.5_r8 * (dp_ptb(i,j,k) / dmg(i,j,k) + dp_ptb(i+1,j,k) / dmg(i+1,j,k)) * &
                 (gz(i+1,j,k) - gz(i,j,k)) / mesh%de_lon(j)
          ! FIXME: Add moisture impact factor.
          du(i,j,k) = du(i,j,k) - (tmp1 + tmp2 + tmp3 + tmp4)
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp1 = 0.5_r8 * (alp_ptb(i,j,k) + alp_ptb(i,j+1,k)) * dmgdy(i,j,k)
          tmp2 = 0.5_r8 * (alp_ref(i,j,k) + alp_ref(i,j+1,k)) * &
                 (p_ptb(i,j+1,k) - p_ptb(i,j,k)) / mesh%de_lat(j)
          tmp3 = (gz_ptb(i,j+1,k) - gz_ptb(i,j,k)) / mesh%de_lat(j)
          tmp4 = 0.5_r8 * (dp_ptb(i,j,k) / dmg(i,j,k) + dp_ptb(i,j+1,k) / dmg(i,j+1,k)) * &
                 (gz(i,j+1,k) - gz(i,j,k)) / mesh%de_lat(j)
          ! FIXME: Add moisture impact factor.
          dv(i,j,k) = dv(i,j,k) - (tmp1 + tmp2 + tmp3 + tmp4)
        end do
      end do
    end do
    end associate

  end subroutine pgf_ptb_run

  subroutine ref_profile_init(this, mesh)

    class(ref_profile_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%mg    , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz    , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz_lev, full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%alp   , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mpt   , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmgdx , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmgdy , full_lon=.true., half_lat=.true., full_lev=.true.)

  end subroutine ref_profile_init

  subroutine ref_profile_clear(this)

    class(ref_profile_type), intent(inout) :: this

    this%mesh => null()
    if (allocated(this%mg    )) deallocate(this%mg    )
    if (allocated(this%gz    )) deallocate(this%gz    )
    if (allocated(this%gz_lev)) deallocate(this%gz_lev)
    if (allocated(this%alp   )) deallocate(this%alp   )
    if (allocated(this%mpt   )) deallocate(this%mpt   )
    if (allocated(this%dmgdx )) deallocate(this%dmgdx )
    if (allocated(this%dmgdy )) deallocate(this%dmgdy )

  end subroutine ref_profile_clear

  subroutine ref_profile_final(this)

    type(ref_profile_type), intent(inout) :: this

    call this%clear()

  end subroutine ref_profile_final

end module pgf_ptb_mod
