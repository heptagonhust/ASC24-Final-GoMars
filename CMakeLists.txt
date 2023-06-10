cmake_minimum_required(VERSION 3.0)

project(gmcore LANGUAGES Fortran C)

option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)
option(ENSURE_ORDER "Ensure computation order in parallel" OFF)
option(USE_BULK_AERO "Use bulk aerosol model" ON)
option(USE_MODAL_AERO "Use modal aerosol model" OFF)
option(USE_CHEM "Use chemistry" OFF)

if (R4)
  add_definitions(-DREAL_KIND=4)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
else ()
  add_definitions(-DREAL_KIND=8)
endif ()
if (ENSURE_ORDER)
  add_definitions(-DENSURE_ORDER)
  set(CPPFLAGS ${CPPFLAGS} -DENSURE_ORDER)
  message(STATUS "Ensure computation order in parallel.")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -fno-range-check -fallow-argument-mismatch -fallow-invalid-boz")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -fcheck=all -g -fbacktrace")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
  endif ()
  # FIXME: OpenMP causes sigmentation faults.
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback -fp-model precise -qopenmp")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O3 -align array64byte -qopenmp") # Add -fp-model precise when checking serial-parallel and restart consistencies.
  endif ()
  if (CMAKE_Fortran_COMPILER_VERSION GREATER 2021.0)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl=sequential")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential")
  endif ()
endif ()

find_package(MPI)

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()
set(EXTERNAL_LIBS netcdff netcdf)

if (DEFINED ENV{GPTL} AND (NOT DEFINED ENV{GPTL_ROOT}))
  set(ENV{GPTL_ROOT} $ENV{GPTL})
endif ()
if (DEFINED ENV{GPTL_ROOT})
  message(STATUS "Found GPTL in $ENV{GPTL_ROOT}.")
  add_definitions(-DHAS_GPTL)
  set(HAS_GPTL TRUE)
  include_directories("$ENV{GPTL_ROOT}/include")
  link_directories("$ENV{GPTL_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS gptlf gptl)
endif ()

if (DEFINED ENV{PIO} AND (NOT DEFINED ENV{PIO_ROOT}))
  set(ENV{PIO_ROOT} $ENV{PIO})
endif ()
if (DEFINED ENV{PIO_ROOT})
  message(STATUS "Found PIO in $ENV{PIO_ROOT}.")
  add_definitions(-DHAS_PIO)
  set(HAS_PIO TRUE)
  include_directories("$ENV{PIO_ROOT}/include")
  link_directories("$ENV{PIO_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS piof pioc)
endif ()

if (DEFINED ENV{LAPACK} AND (NOT DEFINED ENV{LAPACK_ROOT}))
  set(ENV{LAPACK_ROOT} $ENV{LAPACK})
endif ()
if (DEFINED ENV{LAPACK_ROOT})
  message(STATUS "Found LAPACK in $ENV{LAPACK_ROOT}.")
  add_definitions(-DHAS_LAPACK)
  set(HAS_LAPACK TRUE)
  include_directories("$ENV{LAPACK_ROOT}/include")
  link_directories("$ENV{LAPACK_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS lapack)
endif ()

if (DEFINED ENV{ESMF} AND (NOT DEFINED ENV{ESMF_ROOT}))
  set(ENV{ESMF_ROOT} $ENV{ESMF})
endif ()
if (HAS_LAPACK AND DEFINED ENV{ESMF_ROOT})
  message(STATUS "Found ESMF in $ENV{ESMF_ROOT}.")
  add_definitions(-DHAS_ESMF)
  set(HAS_ESMF TRUE)
  include_directories("$ENV{ESMF_ROOT}/include")
  link_directories("$ENV{ESMF_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS esmf)
endif ()

if (DEFINED ENV{MCT} AND (NOT DEFINED ENV{MCT_ROOT}))
  set(ENV{MCT_ROOT} $ENV{MCT})
endif ()
if (DEFINED ENV{MCT_ROOT})
  message(STATUS "Found MCT in $ENV{MCT_ROOT}.")
  add_definitions(-DHAS_MCT)
  set(HAS_MCT TRUE)
  include_directories("$ENV{MCT_ROOT}/include")
  link_directories("$ENV{MCT_ROOT}/lib")
  list(APPEND EXTERNAL_LIBS mct mpeu)
endif ()

if (DEFINED ENV{ECCODES} AND (NOT DEFINED ENV{ECCODES_ROOT}))
  set(ENV{ECCODES_ROOT} $ENV{ECCODES})
endif ()
if (DEFINED ENV{ECCODES_ROOT})
  message(STATUS "Found EcCodes in $ENV{ECCODES_ROOT}.")
  add_definitions(-DHAS_ECCODES)
  include_directories("$ENV{ECCODES_ROOT}/include")
  link_directories("$ENV{ECCODES_ROOT}/lib")
  link_directories("$ENV{ECCODES_ROOT}/lib64")
  list(APPEND EXTERNAL_LIBS eccodes_f90)
else ()
  find_package(ECCODES QUIET)
  if (ECCODES_FOUND)
    list(APPEND EXTERNAL_LIBS eccodes_f90)
  else ()
    message(WARNING "Unable to find EcCodes!")
  endif ()
endif ()

if (DEFINED ENV{MKLROOT})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKLROOT}/include")
  link_directories("$ENV{MKLROOT}/lib")
elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKL_INC}")
  link_directories("$ENV{MKL_LIB}")
endif ()

if (HAS_MKL)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    list(APPEND EXTERNAL_LIBS lapack64)
  endif ()
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

set(gmcore_files
  lib/quadpack.f90
  src/dynamics/adv/adv_batch_mod.F90
  src/dynamics/adv/adv_mod.F90
  src/dynamics/adv/ffsl_mod.F90
  src/dynamics/adv/limiter_mod.F90
  src/dynamics/adv/ppm_mod.F90
  src/dynamics/adv/tvd_mod.F90
  src/dynamics/adv/upwind_mod.F90
  src/dynamics/adv/weno_mod.F90
  src/dynamics/damp/damp_mod.F90
  src/dynamics/damp/div_damp_mod.F90
  src/dynamics/damp/laplace_damp_mod.F90
  src/dynamics/damp/pole_damp_mod.F90
  src/dynamics/damp/smag_damp_mod.F90
  src/dynamics/dynamics_types_mod.F90
  src/dynamics/filter_mod.F90
  src/dynamics/filter_types_mod.F90
  src/dynamics/gmcore_mod.F90
  src/dynamics/interp_mod.F90
  src/dynamics/nh_mod.F90
  src/dynamics/operators_mod.F90
  src/dynamics/pgf/pgf_lin97_mod.F90
  src/dynamics/pgf/pgf_mod.F90
  src/dynamics/pgf/pgf_swm_mod.F90
  src/dynamics/time_schemes_mod.F90
  src/dynamics/tracer_mod.F90
  src/dynamics/tracer_types_mod.F90
  src/dynamics/vert_coord/hybrid_coord_ecmwf_mod.F90
  src/dynamics/vert_coord/hybrid_coord_mod.F90
  src/dynamics/vert_coord/hybrid_coord_ncep_mod.F90
  src/dynamics/vert_coord/hybrid_coord_test_mod.F90
  src/dynamics/vert_coord/mars_vert_coord_mod.F90
  src/dynamics/vert_coord/sigma_coord_mod.F90
  src/dynamics/vert_coord/smooth_coord_mod.F90
  src/dynamics/vert_coord/vert_coord_mod.F90
  src/dynamics/vert_coord/wrf_vert_coord_mod.F90
  src/meshes/cubed_sphere/cubed_sphere_panel_mod.F90
  src/meshes/cubed_sphere/cubed_sphere_mesh_mod.F90
  src/meshes/cubed_sphere/cubed_sphere_array_mod.F90
  src/meshes/latlon/latlon_mesh_mod.F90
  src/meshes/latlon/latlon_halo_mod.F90
  src/meshes/latlon/latlon_parallel_types_mod.F90
  src/meshes/latlon/latlon_parallel_zonal_mod.F90
  src/meshes/latlon/latlon_parallel_mod.F90
  src/meshes/latlon/latlon_decomp_mod.F90
  src/meshes/latlon/latlon_interp_mod.F90
  src/physics/albedo_mod.F90
  src/physics/dp_coupling_mod.F90
  src/physics/mars/mars_cnv_nasa_mod.F90
  src/physics/mars/mars_co2cyc_nasa_mod.F90
  src/physics/mars/mars_data_mod.F90
  src/physics/mars/mars_lsm_nasa_mod.F90
  src/physics/mars/mars_mp_nasa_mod.F90
  src/physics/mars/mars_nasa_mod.F90
  src/physics/mars/mars_pbl_nasa_mod.F90
  src/physics/mars/mars_rad_nasa_mod.F90
  src/physics/mars/mars_rad_rrtmgp_mod.F90
  src/physics/mars/mars_sfc_nasa_mod.F90
  src/physics/pbl/pbl_driver_mod.F90
  src/physics/pbl/pbl_ysu_mod.F90
  src/physics/physics_mod.F90
  src/physics/physics_types_mod.F90
  src/physics/simple_physics_v6.f90
  src/physics/test_forcing_mod.F90
  src/prepare/latlon/latlon_bkg_mod.F90
  src/prepare/latlon/latlon_topo_mod.F90
  src/prepare/prepare_mod.F90
  src/prepare/readers/era5_reader_mod.F90
  src/prepare/readers/cam_reader_mod.F90
  src/prepare/readers/fnl_reader_mod.F90
  src/prepare/readers/mpas_reader_mod.F90
  src/prepare/readers/openmars_reader_mod.F90
  src/prepare/readers/topo_reader_mod.F90
  src/prepare/readers/waccm_reader_mod.F90
  src/prepare/ref_mod.F90
  src/prepare/vert_interp_mod.F90
  src/tests/aquaplanet_test_mod.F90
  src/tests/adv/dcmip12_test_mod.F90
  src/tests/adv/deform_test_mod.F90
  src/tests/adv/moving_vortices_test_mod.F90
  src/tests/adv/solid_rotation_test_mod.F90
  src/tests/baroclinic_wave_test_mod.F90
  src/tests/dcmip31_test_mod.F90
  src/tests/held_suarez_test_mod.F90
  src/tests/ksp15_test_mod.F90
  src/tests/mars_cold_run_mod.F90
  src/tests/mountain_wave_test_mod.F90
  src/tests/rossby_haurwitz_wave_3d_test_mod.F90
  src/tests/steady_state_pgf_test_mod.F90
  src/tests/steady_state_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/shallow_water_waves_test_mod.F90
  src/tests/swm/splash_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  src/tests/swm/vortex_erosion_test_mod.F90
  src/tests/tropical_cyclone_test_mod.F90
  src/utils/accum_mod.F90
  src/utils/allocator_mod.F90
  src/utils/block_mod.F90
  src/utils/const_mod.F90
  src/utils/debug_mod.F90
  src/utils/diag_state_mod.F90
  src/utils/formula_mod.F90
  src/utils/gas_mod.F90
  src/utils/history_mod.F90
  src/utils/initial_mod.F90
  src/utils/math_mod.F90
  src/utils/namelist_mod.F90
  src/utils/perf_mod.F90
  src/utils/process_mod.F90
  src/utils/restart_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/time_mod.F90
  src/utils/time_mod.F90
  src/utils/tridiag_mkl_mod.F90
  src/utils/tridiag_mod.F90
  src/utils/tridiag_spk_mod.F90
)

if (HAS_ESMF)
  list(APPEND gmcore_files
    src/coupler/atm_comp_mod.F90
    src/coupler/lnd_comp_mod.F90
    src/coupler/ocn_comp_mod.F90
    src/coupler/ice_comp_mod.F90
    src/coupler/cpl_comp_mod.F90
  )
endif ()

if (EXISTS ${PROJECT_SOURCE_DIR}/src/physics/ccpp/host/CCPP_CAPS.cmake)
  include(src/physics/ccpp/host/CCPP_CAPS.cmake)
  include(src/physics/ccpp/host/CCPP_SCHEMES.cmake)
  include(src/physics/ccpp/host/CCPP_STATIC_API.cmake)
  list(APPEND gmcore_files ${CAPS} ${SCHEMES} ${API} src/physics/ccpp/host/ccpp_driver_mod.F90)
  if (NOT DEFINED ENV{NCEPLIBS_ROOT})
    message(FATAL_ERROR "NCEPLIBS_ROOT environment variable was not defined!")
  endif ()
  include_directories($ENV{NCEPLIBS_ROOT}/bacio-2.4.1/include_4)
  link_directories($ENV{NCEPLIBS_ROOT}/bacio-2.4.1/lib)
  list(APPEND EXTERNAL_LIBS bacio_4)
  include_directories($ENV{NCEPLIBS_ROOT}/nemsio-2.5.2/include)
  link_directories($ENV{NCEPLIBS_ROOT}/nemsio-2.5.2/lib)
  list(APPEND EXTERNAL_LIBS nemsio)
  include_directories($ENV{NCEPLIBS_ROOT}/w3emc-2.7.3/include_4)
  link_directories($ENV{NCEPLIBS_ROOT}/w3emc-2.7.3/lib)
  list(APPEND EXTERNAL_LIBS w3emc_4)
  include_directories($ENV{NCEPLIBS_ROOT}/w3nco-2.4.1/include_4)
  link_directories($ENV{NCEPLIBS_ROOT}/w3nco-2.4.1/lib)
  list(APPEND EXTERNAL_LIBS w3nco_4)
  include_directories($ENV{NCEPLIBS_ROOT}/g2-3.4.5/include_4)
  link_directories($ENV{NCEPLIBS_ROOT}/g2-3.4.5/lib)
  link_directories($ENV{NCEPLIBS_ROOT}/g2-3.4.5/lib64)
  list(APPEND EXTERNAL_LIBS g2_4)
  include_directories($ENV{NCEPLIBS_ROOT}/sp-2.3.3/include_4)
  link_directories($ENV{NCEPLIBS_ROOT}/sp-2.3.3/lib)
  list(APPEND EXTERNAL_LIBS sp_4)
  add_definitions(-DHAS_CCPP)
endif ()

add_library(gmcore ${gmcore_files})
target_link_libraries(gmcore fortran_container fortran_datetime fiona flogger ${EXTERNAL_LIBS})
if (HAS_GPTL AND HAS_PIO AND HAS_LAPACK AND HAS_ESMF AND HAS_MCT)
  message(STATUS "GPTL, PIO, LAPACK, ESMF, and MCT are found. Building CAM physics and chemistry.")
  include_directories(src/physics/cam/include)
  add_definitions(-DSPMD)
  add_definitions(-DUSE_CSM_SHARE)
  add_definitions(-DPCOLS=16)
  add_definitions(-DPSUBCOLS=1)
  add_definitions(-DPLEV=32)
  add_definitions(-DPCNST=33)
  add_definitions(-DN_RAD_CNST=15)
  add_definitions(-DHAVE_IEEE_ARITHMETIC)
  add_definitions(-DUSE_CONTIGUOUS=)
  add_definitions(-DCLUBB_SGS)
  add_definitions(-DCLUBB_CAM)
  add_definitions(-DCLUBB_REAL_TYPE=8)
  set(USE_CLUBB_SGS ON)
  add_definitions(-DMODAL_AERO_3MODE)
  add_subdirectory(src/physics/cam)
  add_subdirectory(src/chemistry)
  add_definitions(-DHAS_CAM)
else ()
  message(WARNING "GPTL, PIO, LAPACK, ESMF, and/or MCT are not found. Skipping CAM physics and chemistry.")
endif ()

add_executable(gmcore_prepare.exe src/drivers/gmcore_prepare.F90)
target_link_libraries(gmcore_prepare.exe gmcore)

add_executable(gmcore_swm_driver.exe src/drivers/gmcore_swm_driver.F90)
target_link_libraries(gmcore_swm_driver.exe gmcore)

add_executable(gmcore_driver.exe src/drivers/gmcore_driver.F90)
target_link_libraries(gmcore_driver.exe gmcore)

add_executable(gmcore_adv_driver.exe src/drivers/gmcore_adv_driver.F90)
target_link_libraries(gmcore_adv_driver.exe gmcore)
