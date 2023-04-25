! This module connects CCPP to GMCORE.

module ccpp_driver_mod

  use ccpp_static_api
  use CCPP_data, only: GFS_control, GFS_data, GFS_interstitial
  use const_mod
  use namelist_mod, only: hydrostatic, restart, dt_dyn, dt_phys
  use block_mod
  use process_mod
  use time_mod
  use tracer_mod
  use vert_coord_mod

  implicit none

  private

  public ccpp_driver_init
  public ccpp_driver_final
  public ccpp_driver_input_dynamics
  public ccpp_driver_output_physics

contains

  subroutine ccpp_driver_init(namelist_file)

    character(*), intent(in) :: namelist_file

    integer iblk, ithrd, icol
    integer, allocatable :: ncol(:)
    integer idat(8), jdat(8)
    character(:), allocatable, target :: input_nml_file(:)

    call ccpp_driver_final()

    allocate(input_nml_file(1), mold=namelist_file)
    input_nml_file(1) = namelist_file

    allocate(GFS_data(size(blocks)))
    allocate(GFS_interstitial(1)) ! Size 1 for now, will be increased later according to OpenMP thread count.

    idat(:) = 0
    idat(1) = start_time%year
    idat(2) = start_time%month
    idat(3) = start_time%day
    idat(4) = start_time%hour
    idat(5) = start_time%minute
    idat(6) = start_time%second

    jdat(:) = 0
    jdat(1) = curr_time%year
    jdat(2) = curr_time%month
    jdat(3) = curr_time%day
    jdat(4) = curr_time%hour
    jdat(5) = curr_time%minute
    jdat(6) = curr_time%second

    allocate(ncol(size(blocks)))
    do iblk = 1, size(blocks)
      ncol(iblk) = blocks(iblk)%pstate%ncol
    end do

    call GFS_control%init(            &
      nlunit=10                     , &
      fn_nml=namelist_file          , &
      me=proc%id                    , &
      master=0                      , &
      logunit=6                     , &
      isc=1                         , &
      jsc=1                         , &
      nx=sum(ncol)                  , &
      ny=1                          , &
      levs=global_mesh%full_nlev    , &
      cnx=sum(ncol)                 , &
      cny=1                         , &
      gnx=sum(ncol)                 , &
      gny=1                         , &
      dt_dycore=dt_dyn              , &
      dt_phys=dt_phys               , &
      iau_offset=0                  , &
      idat=idat                     , &
      jdat=jdat                     , &
      nwat=ntracers_water           , &
      tracer_names=tracer_names     , &
      tracer_types=tracer_types     , &
      input_nml_file=input_nml_file , &
      tile_num=0                    , &
      blksz=ncol                    , &
      ak=hyai                       , &
      bk=hybi                       , &
      restart=restart               , &
      hydrostatic=hydrostatic       , &
      communicator=proc%comm        , &
      ntasks=proc%np                , &
      nthreads=1                      &
    )

    deallocate(input_nml_file)

    do iblk = 1, size(blocks)
      call GFS_data(iblk)%statein %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%stateout%create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%sfcprop %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%coupling%create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%grid    %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%tbd     %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%cldprop %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%Radtend %create(ncol(iblk), GFS_control)
      call GFS_data(iblk)%intdiag %create(ncol(iblk), GFS_control)
    end do

    do ithrd = 1, size(GFS_interstitial)
      call GFS_interstitial(ithrd)%create(maxval(ncol), GFS_control)
    end do

    do iblk = 1, size(blocks)
      do icol = 1, blocks(iblk)%pstate%ncol
        GFS_data(iblk)%grid%xlon  (icol) = blocks(iblk)%pstate%lon(icol)
        GFS_data(iblk)%grid%xlat  (icol) = blocks(iblk)%pstate%lat(icol)
        GFS_data(iblk)%grid%xlon_d(icol) = blocks(iblk)%pstate%lon(icol) * deg
        GFS_data(iblk)%grid%xlat_d(icol) = blocks(iblk)%pstate%lat(icol) * deg
        GFS_data(iblk)%grid%sinlat(icol) = sin(GFS_data(iblk)%grid%xlat(icol))
        GFS_data(iblk)%grid%coslat(icol) = cos(GFS_data(iblk)%grid%xlat(icol))
        GFS_data(iblk)%grid%area  (icol) = blocks(iblk)%pstate%area(icol)
        GFS_data(iblk)%grid%dx    (icol) = sqrt(GFS_data(iblk)%grid%area(icol))
      end do
    end do

    deallocate(ncol)

  end subroutine ccpp_driver_init

  subroutine ccpp_driver_final()

    if (allocated(GFS_data        )) deallocate(GFS_data        )
    if (allocated(GFS_interstitial)) deallocate(GFS_interstitial)

  end subroutine ccpp_driver_final

  subroutine ccpp_driver_input_dynamics()

    integer iblk, icol, ilev

    do iblk = 1, size(blocks)
      ! Full levels or layers
      do ilev = 1, blocks(iblk)%pstate%nlev
        do icol = 1, blocks(iblk)%pstate%ncol
          GFS_data(iblk)%Statein%ugrs (icol,ilev  ) = blocks(iblk)%pstate%u  (icol,ilev)
          GFS_data(iblk)%Statein%vgrs (icol,ilev  ) = blocks(iblk)%pstate%v  (icol,ilev)
          GFS_data(iblk)%Statein%tgrs (icol,ilev  ) = blocks(iblk)%pstate%t  (icol,ilev)
          GFS_data(iblk)%Statein%vvl  (icol,ilev  ) = blocks(iblk)%pstate%omg(icol,ilev)
          GFS_data(iblk)%Statein%prsl (icol,ilev  ) = blocks(iblk)%pstate%p  (icol,ilev)
          GFS_data(iblk)%Statein%prslk(icol,ilev  ) = blocks(iblk)%pstate%pk (icol,ilev)
          GFS_data(iblk)%Statein%phil (icol,ilev  ) = blocks(iblk)%pstate%z  (icol,ilev) * g
          GFS_data(iblk)%Statein%qgrs (icol,ilev,:) = blocks(iblk)%pstate%q  (icol,ilev,:)
        end do
      end do
      ! Half levels or interfaces
      do ilev = 1, blocks(iblk)%pstate%nlev + 1
        do icol = 1, blocks(iblk)%pstate%ncol
          GFS_data(iblk)%Statein%prsi (icol,ilev) = blocks(iblk)%pstate%p_lev (icol,ilev)
          GFS_data(iblk)%Statein%prsik(icol,ilev) = blocks(iblk)%pstate%pk_lev(icol,ilev)
          GFS_data(iblk)%Statein%phii (icol,ilev) = blocks(iblk)%pstate%z_lev (icol,ilev) * g
        end do
      end do
    end do

  end subroutine ccpp_driver_input_dynamics

  subroutine ccpp_driver_output_physics()

  end subroutine ccpp_driver_output_physics

end module ccpp_driver_mod