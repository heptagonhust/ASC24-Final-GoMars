module bkg_mod

  use flogger
  use namelist_mod
  use const_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use process_mod
  use parallel_mod
  use era5_reader_mod
  use mpas_reader_mod
  use waccm_reader_mod
  use openmars_reader_mod
  use latlon_interp_mod
  use vert_interp_mod
  use tracer_mod
  use operators_mod

  implicit none

  private

  public bkg_read
  public bkg_final
  public bkg_regrid_mgs
  public bkg_calc_mg
  public bkg_calc_ph
  public bkg_regrid_qv
  public bkg_regrid_pt
  public bkg_regrid_u
  public bkg_regrid_v

contains

  subroutine bkg_read(bkg_type_, bkg_file)

    character(*), intent(in) :: bkg_type_
    character(*), intent(in) :: bkg_file

    bkg_type = bkg_type_

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file)
    case ('mpas')
      call mpas_reader_run(bkg_file)
    case ('waccm')
      call waccm_reader_run(bkg_file)
    case ('openmars')
      call openmars_reader_run(bkg_file)
    case default
      if (proc%is_root()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine bkg_read

  subroutine bkg_final()

    call era5_reader_final()
    call mpas_reader_final()
    call waccm_reader_final()
    call openmars_reader_final()

  end subroutine bkg_final

  subroutine bkg_regrid_mgs()

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
      case ('mpas')
        call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_ps, mesh, mgs)
        do_hydrostatic_correct = .false.
      case ('waccm')
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_ps, mesh, p0)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_zs, mesh, z0)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,waccm_nlev), mesh, t0_p)
        call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_t(:,:,waccm_nlev), mesh, t0  )
        do_hydrostatic_correct = .true.
      case ('openmars')
        call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_ps, mesh, mgs)
        do_hydrostatic_correct = .false.
      end select
      ! According to pressure-height formula based on hydrostatic assumption.
      if (do_hydrostatic_correct) then
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            t0(i,j) = t0(i,j) * (p0(i,j) / t0_p(i,j))**lapse_kappa
            mgs(i,j) = p0(i,j) * (1.0_r8 - lapse_rate * (gzs(i,j) / g - z0(i,j)) / t0(i,j))**(1.0_r8 / lapse_kappa)
          end do
        end do
      end if
      call fill_halo(block%halo, mgs, full_lon=.true., full_lat=.true.)
      deallocate(p0, t0, z0, t0_p)
      end associate
    end do

  end subroutine bkg_regrid_mgs

  subroutine bkg_calc_mg()

    integer iblk

    if (proc%is_root()) call log_notice('Calculate dry-air weight.')

    do iblk = 1, size(blocks)
      call calc_mg (blocks(iblk), blocks(iblk)%dstate(1))
      call calc_dmg(blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine bkg_calc_mg

  subroutine bkg_calc_ph()

    integer iblk

    if (proc%is_root()) call log_notice('Calculating full pressure.')

    do iblk = 1, size(blocks)
      call calc_ph(blocks(iblk), blocks(iblk)%dstate(1))
    end do

  end subroutine bkg_calc_ph

  subroutine bkg_regrid_pt()

    real(r8), pointer :: qv(:,:,:)
    real(r8), allocatable, dimension(:,:,:) :: t1, pt1, p1
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid temperature and calculate modified potential temperature.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 mg    => blocks(iblk)%dstate(1)%mg, & ! in
                 ph    => blocks(iblk)%dstate(1)%ph, & ! in
                 t     => blocks(iblk)%dstate(1)%t , & ! out
                 tv    => blocks(iblk)%dstate(1)%tv, & ! out
                 pt    => blocks(iblk)%dstate(1)%pt)   ! out
        call tracer_get_array(iblk, idx_qv, qv)
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
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg(i,j,1:mesh%full_nlev), t(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv(i,j,:) = virtual_temperature(t(i,j,:), qv(i,j,:), qv(i,j,:))
              pt(i,j,:) = modified_potential_temperature(t(i,j,:), ph(i,j,:), qv(i,j,:))
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
              call vert_interp_log_linear(p1(i,j,:), pt1(i,j,:), mg(i,j,1:mesh%full_nlev), pt(i,j,1:mesh%full_nlev), allow_extrap=.false.)
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
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg(i,j,1:mesh%full_nlev), t(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv(i,j,:) = virtual_temperature(t(i,j,:), 0.0_r8, 0.0_r8)
              pt(i,j,:) = modified_potential_temperature(t(i,j,:), ph(i,j,:), 0.0_r8)
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
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), mg(i,j,1:mesh%full_nlev), t(i,j,1:mesh%full_nlev), allow_extrap=.true.)
              tv(i,j,:) = virtual_temperature(t(i,j,:), 0.0_r8, 0.0_r8)
              pt(i,j,:) = modified_potential_temperature(t(i,j,:), ph(i,j,:), 0.0_r8)
            end do
          end do
          deallocate(t1, p1)
        end select
        call fill_halo(block%filter_halo, pt, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
      end associate
    end do

  end subroutine bkg_regrid_pt

  subroutine bkg_regrid_u()

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
          call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_u (:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg(i,j,1:mesh%full_nlev), u(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u1, p1)
      case ('mpas')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,mpas_nlev))
        do k = 1, mpas_nlev
          call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg(i,j,1:mesh%full_nlev), u(i,j,1:mesh%full_nlev), allow_extrap=.false.)
          end do
        end do
        deallocate(u1, p1)
      case ('waccm')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,waccm_nlev))
        do k = 1, waccm_nlev
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
              call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg(i,j,1:mesh%full_nlev), u(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u1, p1)
      case ('openmars')
        allocate(u1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        allocate(p1(mesh%half_ims:mesh%half_ime,mesh%full_jms:mesh%full_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%half_ids, mesh%half_ide
            call vert_interp_linear(p1(i,j,:), u1(i,j,:), mg(i,j,1:mesh%full_nlev), u(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(u1, p1)
      end select
      call fill_halo(block%halo, u, full_lon=.false., full_lat=.true., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_u

  subroutine bkg_regrid_v()

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
          call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_v (:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_pd(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg(i,j,1:mesh%full_nlev), v(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v1, p1)
      case ('mpas')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mpas_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,mpas_nlev))
        do k = 1, mpas_nlev
          call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg(i,j,1:mesh%full_nlev), v(i,j,1:mesh%full_nlev), allow_extrap=.false.)
          end do
        end do
        deallocate(v1, p1)
      case ('waccm')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,waccm_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,waccm_nlev))
        do k = 1, waccm_nlev
          call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg(i,j,1:mesh%full_nlev), v(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v1, p1)
      case ('openmars')
        allocate(v1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        allocate(p1(mesh%full_ims:mesh%full_ime,mesh%half_jms:mesh%half_jme,openmars_nlev))
        do k = 1, openmars_nlev
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_v(:,:,k), mesh, v1(:,:,k))
          call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            call vert_interp_linear(p1(i,j,:), v1(i,j,:), mg(i,j,1:mesh%full_nlev), v(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(v1, p1)
      end select
      call fill_halo(block%halo, v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_v

  subroutine bkg_regrid_qv()

    real(r8), pointer :: qv(:,:,:)
    real(r8), allocatable, dimension(:,:,:) :: q1, p1
    integer iblk, i, j, k

    if (proc%is_root()) call log_notice('Regrid water vapor mixing ratio.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%filter_mesh , &
                 mg    => blocks(iblk)%dstate(1)%mg)   ! in
      call tracer_get_array(iblk, idx_qv, qv)
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
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg(i,j,1:mesh%full_nlev), qv(i,j,1:mesh%full_nlev), allow_extrap=.true.)
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
            call vert_interp_linear(p1(i,j,:), q1(i,j,:), mg(i,j,1:mesh%full_nlev), qv(i,j,1:mesh%full_nlev), allow_extrap=.true.)
          end do
        end do
        deallocate(q1, p1)
      end select
      call fill_halo(block%filter_halo, qv, full_lon=.true., full_lat=.true., full_lev=.true., cross_pole=.true.)
      end associate
    end do

  end subroutine bkg_regrid_qv

end module bkg_mod
