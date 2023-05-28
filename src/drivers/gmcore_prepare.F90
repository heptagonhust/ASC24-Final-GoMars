program gmcore_prepare

  use fiona
  use string
  use time_mod
  use initial_mod
  use namelist_mod
  use gmcore_mod
  use prepare_mod

  implicit none

  character(256) namelist_file

  call get_command_argument(1, namelist_file)

  call parse_namelist(namelist_file)

  time_scheme = 'N/A'

  call fiona_init()

  call gmcore_init_stage1(namelist_file)
  call prepare_topo()
  call gmcore_init_stage2(namelist_file)
  call prepare_bkg()

  if (initial_file == 'N/A') then
    write(initial_file, '("gmcore.", A, "x", A, "x", A, ".i0.nc")') &
      to_str(nlon), &
      to_str(nlat), &
      to_str(nlev)
  end if
  call initial_write(initial_file, initial_time)

  call prepare_final()

  call process_final()

end program gmcore_prepare
