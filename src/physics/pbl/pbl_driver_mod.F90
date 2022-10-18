module pbl_driver_mod

  use const_mod
  use namelist_mod
  use physics_types_mod
  use pbl_ysu_mod

  implicit none

  private

  public pbl_init
  public pbl_run
  public pbl_final

contains

  subroutine pbl_init()

  end subroutine pbl_init

  subroutine pbl_run(pstate, ptend, dt)

    type(pstate_type), intent(inout) :: pstate
    type(ptend_type ), intent(inout) :: ptend
    real(8), intent(in) :: dt

    logical flag_qi
    integer ysu_topdown_pblmix

    select case (pbl_scheme)
    case ('ysu')
      call ysu(                                 &
        u3d               =pstate%u           , &
        v3d               =pstate%v           , &
        th3d              =pstate%pt          , &
        t3d               =pstate%t           , &
        qv3d              =pstate%qv          , &
        qc3d              =pstate%qc          , &
        qi3d              =pstate%qi          , &
        p3d               =pstate%p           , &
        p3di              =pstate%p_lev       , &
        pi3d              =pstate%p_exn       , &
        rublten           =ptend%dudt         , &
        rvblten           =ptend%dvdt         , &
        rthblten          =ptend%dptdt        , &
        rqvblten          =ptend%dqvdt        , &
        rqcblten          =ptend%dqcdt        , &
        rqiblten          =ptend%dqidt        , &
        flag_qi           =flag_qi            , &
        cp                =cpd                , &
        g                 =g                  , &
        rovcp             =Rd_o_cpd           , &
        rd                =Rd                 , &
        rovg              =Rd_o_g             , &
        ep1               =Rv_o_Rd - 1        , &
        ep2               =Rd_o_Rv            , &
        karman            =karman             , &
        xlv               =Lv                 , &
        rv                =Rv                 , &
        dz8w              =pstate%dz          , &
        psfc              =pstate%ps          , &
        znt               =pstate%z0          , &
        ust               =pstate%ustar       , &
        hpbl              =pstate%pblh        , &
        psim              =pstate%psim        , &
        psih              =pstate%psih        , &
        xland             =pstate%land        , &
        hfx               =pstate%hfx         , &
        qfx               =pstate%qfx         , &
        wspd              =pstate%wsp         , &
        br                =pstate%Rib         , &
        dt                =dt                 , &
        kpbl2d            =pstate%pblk        , &
        exch_h            =pstate%exch_h      , &
        wstar             =pstate%wstar       , &
        delta             =pstate%delta       , &
        u10               =pstate%u10         , &
        v10               =pstate%v10         , &
        uoce              =pstate%uos         , &
        voce              =pstate%vos         , &
        rthraten          =ptend%dptdt_rad    , &
        ysu_topdown_pblmix=ysu_topdown_pblmix , &
        ids=1, ide=pstate%ncol, jds=1, jde=1, kds=1, kde=pstate%nlev, &
        ims=1, ime=pstate%ncol, jms=1, jme=1, kms=1, kme=pstate%nlev, &
        its=1, ite=pstate%ncol, jts=1, jte=1, kts=1, kte=pstate%nlev  &
        !ctopo             =...                , &
        !ctopo2            =...                , &
      )
    end select

  end subroutine pbl_run

  subroutine pbl_final()

  end subroutine pbl_final

end module pbl_driver_mod
