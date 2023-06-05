module prepare_mod

  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use process_mod
  use block_mod
  use topo_reader_mod
  use latlon_topo_mod
  use latlon_bkg_mod
  use ref_mod
  use operators_mod
  use tracer_mod

  implicit none

contains

  subroutine prepare_topo()

    integer iblk, i, j

    if (.not. use_aqua_planet) then
      call topo_reader_run(topo_file, min_lon, max_lon, min_lat, max_lat)
      do iblk = 1, size(blocks)
        call latlon_topo_regrid(blocks(iblk))
      end do
      if (use_topo_smooth) then
        do iblk = 1, size(blocks)
          call latlon_topo_smooth(blocks(iblk))
        end do
      end if
      call ref_calc_ps()
    end if

    do iblk = 1, size(blocks)
      associate (mesh    => blocks(iblk)%mesh          , &
                 gzs     => blocks(iblk)%static%gzs    , & ! in
                 dzsdlon => blocks(iblk)%static%dzsdlon, & ! out
                 dzsdlat => blocks(iblk)%static%dzsdlat)   ! out
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          dzsdlon(i,j) = (gzs(i+1,j) - gzs(i,j)) / g / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          dzsdlat(i,j) = (gzs(i,j+1) - gzs(i,j)) / g / mesh%de_lat(j)
        end do
      end do
      call fill_halo(blocks(iblk)%halo, dzsdlon, full_lon=.false., full_lat=.true.)
      call fill_halo(blocks(iblk)%halo, dzsdlat, full_lon=.true., full_lat=.false.)
      end associate
    end do

  end subroutine prepare_topo

  subroutine prepare_bkg()

    integer iblk

    call latlon_bkg_read(min_lon, max_lon, min_lat, max_lat)

    call latlon_bkg_regrid_mgs()
    call latlon_bkg_calc_mg()
    call latlon_bkg_regrid_qv()
    call latlon_bkg_calc_ph()
    call latlon_bkg_regrid_pt()
    call latlon_bkg_regrid_u()
    call latlon_bkg_regrid_v()

    do iblk = 1, size(blocks)
      call calc_gz_lev(blocks(iblk), blocks(iblk)%dstate(1))
      call tracer_calc_qm(blocks(iblk))
    end do

  end subroutine prepare_bkg

  subroutine prepare_final()

    call topo_reader_final()
    call latlon_bkg_final()

  end subroutine prepare_final

end module prepare_mod
