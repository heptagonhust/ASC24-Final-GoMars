module pbl_driver_mod

  use const_mod
  use namelist_mod
  use tracer_types_mod
  use wrf_physics_types_mod
  use pbl_ysu_mod

  implicit none

  private

  public pbl_init
  public pbl_run
  public pbl_final

contains

  subroutine pbl_init()

  end subroutine pbl_init

  subroutine pbl_run(state, tend, dt)

    type(wrf_state_type), intent(inout) :: state
    type(wrf_tend_type ), intent(inout) :: tend
    real(r8), intent(in) :: dt

    logical flag_qi
    integer ysu_topdown_pblmix

    select case (pbl_scheme)
    case ('ysu')
      call ysu(                                   &
        u3d               =state%u              , & ! done
        v3d               =state%v              , & ! done
        th3d              =state%pt             , & ! no use
        t3d               =state%t              , & ! done
        qv3d              =state%q(:,:,idx_qv)  , & ! done
        qc3d              =state%q(:,:,idx_qc)  , & ! no data
        qi3d              =state%q(:,:,idx_qi)  , & ! no data
        p3d               =state%p              , & ! done
        p3di              =state%p_lev          , & ! done
        pi3d              =state%pk             , & ! done
        rublten           =tend%dudt            , & ! out
        rvblten           =tend%dvdt            , & ! out
        rthblten          =tend%dptdt           , & ! out
        rqvblten          =tend%dqdt(:,:,idx_qv), & ! out
        rqcblten          =tend%dqdt(:,:,idx_qc), & ! out
        rqiblten          =tend%dqdt(:,:,idx_qi), & ! out
        flag_qi           =flag_qi              , & ! no use
        cp                =cpd                  , & ! done
        g                 =g                    , & ! done
        rovcp             =rd_o_cpd             , & ! done
        rd                =rd                   , & ! done
        rovg              =rd_o_g               , & ! done
        ep1               =rv_o_rd - 1          , & ! done
        ep2               =rd_o_rv              , & ! done
        karman            =ka                   , & ! done
        xlv               =lv                   , & ! done
        rv                =rv                   , & ! done
        dz8w              =state%dz             , & ! done
        psfc              =state%ps             , & ! done
        znt               =state%z0             , & ! ?
        ust               =state%ustar          , & ! ?
        psim              =state%psim           , & ! ?
        psih              =state%psih           , & ! ?
        xland             =state%land           , & ! done
        hfx               =state%hfx            , & ! ?
        qfx               =state%qfx            , & ! ?
        wspd              =state%wsp_bot        , & ! done
        br                =state%br             , & ! ?
        dt                =dt                   , & ! done
        hpbl              =state%pblh           , & ! out
        kpbl2d            =state%pblk           , & ! out
        exch_h            =state%exch_h         , & ! out
        wstar             =state%wstar          , & ! out
        delta             =state%delta          , & ! out
        u10               =state%u10            , & ! ?
        v10               =state%v10            , & ! ?
        uoce              =state%uos            , & ! ?
        voce              =state%vos            , & ! ?
        rthraten          =tend%dptdt_rad       , & ! ?
        ysu_topdown_pblmix=ysu_topdown_pblmix   , &
        ids=1, ide=state%mesh%ncol, jds=1, jde=1, kds=1, kde=state%mesh%nlev, &
        ims=1, ime=state%mesh%ncol, jms=1, jme=1, kms=1, kme=state%mesh%nlev, &
        its=1, ite=state%mesh%ncol, jts=1, jte=1, kts=1, kte=state%mesh%nlev  &
        !ctopo             =...                , &
        !ctopo2            =...                , &
      )
    end select

  end subroutine pbl_run

  subroutine pbl_final()

  end subroutine pbl_final

end module pbl_driver_mod
