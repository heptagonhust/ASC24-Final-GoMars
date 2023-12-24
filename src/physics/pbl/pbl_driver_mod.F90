module pbl_driver_mod

  use const_mod
  use namelist_mod
  use tracer_types_mod
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
    real(r8), intent(in) :: dt

    logical flag_qi
    integer ysu_topdown_pblmix

    select case (pbl_scheme)
    case ('ysu')
      call ysu(                                    &
        u3d               =pstate%u              , & ! done
        v3d               =pstate%v              , & ! done
        th3d              =pstate%pt             , & ! no use
        t3d               =pstate%t              , & ! done
        qv3d              =pstate%qv             , & ! done
        qc3d              =pstate%qc             , & ! no data
        qi3d              =pstate%qi             , & ! no data
        p3d               =pstate%p              , & ! done
        p3di              =pstate%p_lev          , & ! done
        pi3d              =pstate%pk             , & ! done
        rublten           =ptend%dudt            , & ! out
        rvblten           =ptend%dvdt            , & ! out
        rthblten          =ptend%dptdt           , & ! out
        rqvblten          =ptend%dqdt(:,:,idx_qv), & ! out
        rqcblten          =ptend%dqdt(:,:,idx_qc), & ! out
        rqiblten          =ptend%dqdt(:,:,idx_qi), & ! out
        flag_qi           =flag_qi               , & ! no use
        cp                =cpd                   , & ! done
        g                 =g                     , & ! done
        rovcp             =rd_o_cpd              , & ! done
        rd                =rd                    , & ! done
        rovg              =rd_o_g                , & ! done
        ep1               =rv_o_rd - 1           , & ! done
        ep2               =rd_o_rv               , & ! done
        karman            =ka                    , & ! done
        xlv               =lv                    , & ! done
        rv                =rv                    , & ! done
        dz8w              =pstate%dz             , & ! done
        psfc              =pstate%ps             , & ! done
        znt               =pstate%z0             , & ! ?
        ust               =pstate%ustar          , & ! ?
        psim              =pstate%psim           , & ! ?
        psih              =pstate%psih           , & ! ?
        xland             =pstate%land           , & ! done
        hfx               =pstate%hfx            , & ! ?
        qfx               =pstate%qfx            , & ! ?
        wspd              =pstate%wsb            , & ! done
        br                =pstate%rib            , & ! ?
        dt                =dt                    , & ! done
        hpbl              =pstate%pblh           , & ! out
        kpbl2d            =pstate%pblk           , & ! out
        exch_h            =pstate%exch_h         , & ! out
        wstar             =pstate%wstar          , & ! out
        delta             =pstate%delta          , & ! out
        u10               =pstate%u10            , & ! ?
        v10               =pstate%v10            , & ! ?
        uoce              =pstate%uos            , & ! ?
        voce              =pstate%vos            , & ! ?
        rthraten          =ptend%dptdt_rad       , & ! ?
        ysu_topdown_pblmix=ysu_topdown_pblmix    , &
        ids=1, ide=pstate%mesh%ncol, jds=1, jde=1, kds=1, kde=pstate%mesh%nlev, &
        ims=1, ime=pstate%mesh%ncol, jms=1, jme=1, kms=1, kme=pstate%mesh%nlev, &
        its=1, ite=pstate%mesh%ncol, jts=1, jte=1, kts=1, kte=pstate%mesh%nlev  &
        !ctopo             =...                , &
        !ctopo2            =...                , &
      )
    end select

  end subroutine pbl_run

  subroutine pbl_final()

  end subroutine pbl_final

end module pbl_driver_mod
