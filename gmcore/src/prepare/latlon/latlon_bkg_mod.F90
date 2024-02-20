! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_bkg_mod

  use flogger
  use namelist_mod
  use const_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use process_mod
  use latlon_parallel_mod
  use era5_reader_mod
  use cam_reader_mod
  use mpas_reader_mod
  use waccm_reader_mod
  use openmars_reader_mod
  use latlon_interp_mod
  use vert_interp_mod
  use tracer_mod
  use operators_mod

  implicit none

  private

  public latlon_bkg_read
  public latlon_bkg_final
  public latlon_bkg_regrid_mgs
  public latlon_bkg_calc_mg
  public latlon_bkg_calc_ph
  public latlon_bkg_regrid_qv
  public latlon_bkg_regrid_qc
  public latlon_bkg_regrid_qi
  public latlon_bkg_regrid_nc
  public latlon_bkg_regrid_ni
  public latlon_bkg_regrid_pt
  public latlon_bkg_regrid_u
  public latlon_bkg_regrid_v

contains

  subroutine latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)

    real(r8), intent(in) :: min_lon, max_lon, min_lat, max_lat

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file, min_lon, max_lon, min_lat, max_lat)
    case ('cam')
      call cam_reader_run(bkg_file, min_lon, max_lon, min_lat, max_lat)
    case ('mpas')
      call mpas_reader_run(bkg_file)
    case ('waccm')
      call waccm_reader_run(bkg_file)
    case ('openmars')
      call openmars_reader_run(bkg_file)
    case default
      if (proc%is_root()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine latlon_bkg_read

  subroutine latlon_bkg_final()

    call era5_reader_final()
    call mpas_reader_final()
    call waccm_reader_final()
    call openmars_reader_final()

  end subroutine latlon_bkg_final

  subroutine latlon_bkg_regrid_mgs()

    real(r8), allocatable, dimension(:,:) :: p0, t0, z0, t0_p
    real(r8) lapse_kappa
    integer iblk, i, j

    logical do_hydrostatic_correct

    if (proc%is_root()) call log_notice('Regrid surface dry-air weight.')

    lapse_kappa = lapse_rate * rd_o_g
    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)              , &
                 mesh  => blocks(iblk)%mesh         , &
                 gzs   => blocks(iblk)%static%gzs   , &
                 mgs   => blocks(iblk)%dstate(1)%mgs)
      allocate(p0  (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
      allocate(t0  (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
      allocate(z0  (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))
      allocate(t0_p(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme))

      select case (bkg_type)
      case ('era5')
        call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_psd, mesh, p0)
        call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_zs, mesh, z0)
        call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,era5_nlev), mesh, t0)
        t0_p = era5_lev(era5_nlev)
        do_hydrostatic_correct = .true.
      case ('cam')
        call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_ps, mesh, mgs%d)
        do_hydrostatic_correct = .false.
      case ('mpas')
        call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_ps, mesh, mgs%d)
        do_hydrostatic_correct = .false.
      case ('waccm')
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_ps, mesh, p0)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_zs, mesh, z0)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,waccm_nlev), mesh, t0_p)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_t(:,:,waccm_nlev), mesh, t0  )
        do_hydrostatic_correct = .true.
      case ('openmars')
        call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_ps, mesh, mgs%d)
        do_hydrostatic_correct = .false.
      end select
      ! According to pressure-height formula based on hydrostatic assumption.
      if (do_hydrostatic_correct) then
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            t0(i,j) = t0(i,j) * (p0(i,j) / t0_p(i,j))**lapse_kappa
            mgs%d(i,j) = p0(i,j) * (1.0_r8 - lapse_rate * (gzs%d(i,j) / g - z0(i,j)) / t0(i,j))**(1.0_r8 / lapse_kappa)
          end do
        end do
      end if
      call fill_halo(mgs)
      deallocate(p0, t0, z0, t0_p)
      end associate
    end do

  end subroutine latlon_bkg_regrid_mgs

  subroutine latlon_bkg_calc_mg()

    integer iblk

    if (proc%is_root()) call log_notice('Calculate dry-air weight.')

    do iblk = 1, size(blocks)
      call calc_mg (blocks(iblk), blocks(iblk)%dstate(1))
      call calc_dmg(blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine latlon_bkg_calc_mg

  subroutine latlon_bkg_calc_ph()

    integer iblk

    if (proc%is_root()) call log_notice('Calculating full pressure.')

    do iblk = 1, size(blocks)
      call calc_ph(blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine latlon_bkg_calc_ph

  subroutine latlon_bkg_regrid_pt()

    real(r8), allocatable, dimension(:,:,:) :: t1, pt1, p1
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid temperature and calculate modified potential temperature.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 ph    => blocks(iblk)%dstate(1)%ph, & ! in
                 q     => tracers(iblk)%q          , & ! in
                 qm    => tracers(iblk)%qm         , & ! in
                 t     => blocks(iblk)%dstate(1)%t , & ! out
                 tv    => blocks(iblk)%dstate(1)%tv, & ! out
                 pt    => blocks(iblk)%dstate(1)%pt)   ! out
        select case (bkg_type)
        case ('era5')
          allocate(t1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
          allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
          do k = 1, era5_nlev
            call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t (:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv%d(i,j,:) = virtual_temperature(t%d(i,j,:), q%d(i,j,:,idx_qv), qm%d(i,j,:))
              pt%d(i,j,:) = modified_potential_temperature(t%d(i,j,:), ph%d(i,j,:), q%d(i,j,:,idx_qv))
            end do
          end do
          deallocate(t1, p1)
        case ('cam')
          allocate(t1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
          allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
          do k = 1, cam_nlev
            call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_t(:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv%d(i,j,:) = virtual_temperature(t%d(i,j,:), q%d(i,j,:,idx_qv), qm%d(i,j,:))
              pt%d(i,j,:) = modified_potential_temperature(t%d(i,j,:), ph%d(i,j,:), q%d(i,j,:,idx_qv))
            end do
          end do
          deallocate(t1, p1)
        case ('mpas')
          allocate(pt1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
          allocate(p1 (mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
          do k = 1, mpas_nlev
            call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_pt(:,:,k), mesh, pt1(:,:,k))
            call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_p (:,:,k), mesh, p1 (:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p1(i,j,:), pt1(i,j,:), mg%d(i,j,1:mesh%full_nlev), pt%d(i,j,1:mesh%full_nlev), allow_extrap=.false.)
            end do
          end do
          deallocate(pt1, p1)
        case ('waccm')
          allocate(t1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
          allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
          do k = 1, waccm_nlev
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_t(:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv%d(i,j,:) = virtual_temperature(t%d(i,j,:), 0.0_r8, 0.0_r8)
              pt%d(i,j,:) = modified_potential_temperature(t%d(i,j,:), ph%d(i,j,:), 0.0_r8)
            end do
          end do
          deallocate(t1, p1)
        case ('openmars')
          allocate(t1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
          allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
          do k = 1, openmars_nlev
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_t(:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg%d(i,j,1:mesh%full_nlev), t%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv%d(i,j,:) = virtual_temperature(t%d(i,j,:), 0.0_r8, 0.0_r8)
              pt%d(i,j,:) = modified_potential_temperature(t%d(i,j,:), ph%d(i,j,:), 0.0_r8)
            end do
          end do
          deallocate(t1, p1)
        end select
        call fill_halo(pt)
      end associate
    end do

  end subroutine latlon_bkg_regrid_pt

  subroutine latlon_bkg_regrid_u()

    real(r8), allocatable, dimension(:,:,:) :: u1, p1
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid u wind component.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 mg    => blocks(iblk)%dstate(1)%mg, &
                 u     => blocks(iblk)%dstate(1)%u_lon)
      select case (bkg_type)
      case ('era5')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_u (:,:,k), mesh, u1(:,:,k), extrap=.true., zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k), extrap=.true.)
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('cam')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_lon_edge(cam_lon, cam_slat, cam_us(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(cam_lon, cam_lat , cam_p (:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('mpas')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
        do k = 1, mpas_nlev
          call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.false.)
          end do
        end do
      case ('waccm')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        do k = 1, waccm_nlev
          call latlon_interp_bilinear_lon_edge(waccm_lon, waccm_lat, waccm_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
              call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('openmars')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg%d(i,j,1:mesh%full_nlev), u%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      end select
      deallocate(u1, p1)
      call fill_halo(u)
      end associate
    end do

  end subroutine latlon_bkg_regrid_u

  subroutine latlon_bkg_regrid_v()

    real(r8), allocatable, dimension(:,:,:) :: v1, p1
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid v wind component.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 mg    => blocks(iblk)%dstate(1)%mg, &
                 v     => blocks(iblk)%dstate(1)%v_lat)
      select case (bkg_type)
      case ('era5')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,era5_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_v (:,:,k), mesh, v1(:,:,k), extrap=.true.)
          call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k), extrap=.true.)
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('cam')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_lat_edge(cam_slon, cam_lat, cam_vs(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(cam_lon , cam_lat, cam_p (:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('mpas')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mpas_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mpas_nlev))
        do k = 1, mpas_nlev
          call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.false.)
          end do
        end do
      case ('waccm')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,waccm_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,waccm_nlev))
        do k = 1, waccm_nlev
          call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      case ('openmars')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg%d(i,j,1:mesh%full_nlev), v%d(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
      end select
      deallocate(v1, p1)
      call fill_halo(v)
      end associate
    end do

  end subroutine latlon_bkg_regrid_v

  subroutine latlon_bkg_regrid_qv()

    real(r8), allocatable, dimension(:,:,:) :: q1, p1
    integer iblk, i, j, k

    if (idx_qv == 0) return

    if (proc%is_root()) call log_notice('Regrid water vapor mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 q     => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qv(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qv), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      case ('cam')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_q(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qv), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      case ('waccm')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        do k = 1, waccm_nlev
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_q(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qv), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      end select
      call fill_halo(q, idx_qv)
      end associate
    end do

  end subroutine latlon_bkg_regrid_qv

  subroutine latlon_bkg_regrid_qc()

    real(r8), allocatable, dimension(:,:,:) :: q1, p1
    integer iblk, i, j, k

    if (idx_qc == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud liquid mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 q     => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_ql(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qc), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      case ('cam')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_cldliq(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qc), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      end select
      call fill_halo(q, idx_qc)
      end associate
    end do

  end subroutine latlon_bkg_regrid_qc

  subroutine latlon_bkg_regrid_qi()

    real(r8), allocatable, dimension(:,:,:) :: q1, p1
    integer iblk, i, j, k

    if (idx_qi == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud ice mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 q     => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('era5')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,era5_nlev))
        do k = 1, era5_nlev
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_qi(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qi), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      case ('cam')
        allocate(q1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_cldice(:,:,k), mesh, q1(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_qi), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      end select
      call fill_halo(q, idx_qi)
      end associate
    end do

  end subroutine latlon_bkg_regrid_qi

  subroutine latlon_bkg_regrid_nc()

    real(r8), allocatable, dimension(:,:,:) :: n1, p1
    integer iblk, i, j, k

    if (idx_nc == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud liquid number mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 q     => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('cam')
        allocate(n1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_numliq(:,:,k), mesh, n1(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), n1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_nc), allow_extrap=.true.)
          end do
        end do
        deallocate(n1, p1)
      end select
      call fill_halo(q, idx_nc)
      end associate
    end do

  end subroutine latlon_bkg_regrid_nc

  subroutine latlon_bkg_regrid_ni()

    real(r8), allocatable, dimension(:,:,:) :: n1, p1
    integer iblk, i, j, k

    if (idx_ni == 0) return

    if (proc%is_root()) call log_notice('Regrid cloud ice number mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 q     => tracers(iblk)%q          )   ! out
      select case (bkg_type)
      case ('cam')
        allocate(n1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%full_jms:mesh%full_jme,cam_nlev))
        do k = 1, cam_nlev
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_numice(:,:,k), mesh, n1(:,:,k))
          call latlon_interp_bilinear_cell(cam_lon, cam_lat, cam_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), n1(i,j,:), mg%d(i,j,1:mesh%full_nlev), q%d(i,j,1:mesh%full_nlev,idx_ni), allow_extrap=.true.)
          end do
        end do
        deallocate(n1, p1)
      end select
      call fill_halo(q, idx_ni)
      end associate
    end do

  end subroutine latlon_bkg_regrid_ni

end module latlon_bkg_mod
