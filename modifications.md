# Introduction to optimization methods 

## Compile time optimization
- Automatic vectorization
- IPO (Interprocedural Optimization)
- Static Linking: link all libraries to enhance performance. This change enables more efficient compiler and linker optimizations and reduces the overhead associated with library calls.

## Communication and Parallelization
- Communication-Computation Overlap: Implemented asynchronous communication interface based on MPI, overlapping partial communication and computation to exploit the potential of parallel computing.
- MPI Process Splitting Optimization: Proper splitting in both north and south directions to reduce the maximum communication volume of single MPI_ISEND/RECV, thereby reducing communication latency.
- MPI & OpenMP Hybrid Parallel Optimization: Implemented hybrid parallelization of MPI and OpenMP, enhancing the flexibility of parallel computing and adaptability to different grid resolutions.
- Reduced Redundant Communication: Grid regions of partial communication are not computed or updated, reducing unnecessary communication in those areas.

## Numa System
- Numa system optimization: Utilized specific binding strategies and tools like Slurm, MPI, and numactl to manage processor allocations efficiently, which eliminates memory access delays and reduces cache misses.


 
# Code Modification Details 
**CMakeLists.txt**
```cmake
diff --git a/gmcore/CMakeLists.txt b/gmcore/CMakeLists.txt
index b07b9ae..23b6739 100644
--- a/gmcore/CMakeLists.txt
+++ b/gmcore/CMakeLists.txt
@@ -56,8 +56,7 @@ if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
   endif ()
   # FIXME: OpenMP causes sigmentation faults.
   # list(APPEND fortran_flags -fopenmp)
-elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
-# elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
+elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
   add_definitions(-DFC_IS_INTEL)
   if (CMAKE_BUILD_TYPE STREQUAL "Debug")
     list(APPEND fortran_flags
@@ -73,8 +72,7 @@ elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID
       -no-wrap-margin
       -O3
       -align array64byte
-      # -ipo
-      # -qopenmp
+      -qopenmp
       # -fp-model precise
     )
   endif ()
@@ -86,10 +84,6 @@ elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID
 endif ()
 
 list(APPEND fortran_flags -heap-arrays)
-# list(APPEND fortran_flags -fstack-arrays)
-list(APPEND fortran_flags -xCORE-AVX512)
-
-
 string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${fortran_flags}")
 
 
@@ -115,23 +109,7 @@ else ()
     message(FATAL_ERROR "Unable to find pkg-config library!")
   endif ()
 endif ()
-
-if (DEFINED ENV{H5DIR})
-include_directories("$ENV{H5DIR}/include")
-link_directories("$ENV{H5DIR}/lib")
-endif ()
-if (DEFINED ENV{CURLDIR})
-include_directories("$ENV{CURLDIR}/include")
-link_directories("$ENV{CURLDIR}/lib")
-endif ()
-if (DEFINED ENV{XML2DIR})
-include_directories("$ENV{XML2DIR}/include")
-link_directories("$ENV{XML2DIR}/lib")
-endif ()
-
-set(EXTERNAL_LIBS netcdff -lhdf5_hl -lhdf5 -lm -lz -lsz -lzstd -lxml2 -lcurl)
-set(EXTERNAL_LIBS netcdf -lhdf5_hl -lhdf5 -lm -lz -lsz -lzstd -lxml2 -lcurl)
-#set(EXTERNAL_LIBS netcdff netcdf)
+set(EXTERNAL_LIBS netcdff netcdf)
 
 if (DEFINED ENV{GPTL} AND (NOT DEFINED ENV{GPTL_ROOT}))
   set(ENV{GPTL_ROOT} $ENV{GPTL})
@@ -142,8 +120,7 @@ if (DEFINED ENV{GPTL_ROOT})
   set(HAS_GPTL TRUE)
   include_directories("$ENV{GPTL_ROOT}/include")
   link_directories("$ENV{GPTL_ROOT}/lib")
-  list(APPEND EXTERNAL_LIBS gptlf -lm)
-  list(APPEND EXTERNAL_LIBS gptl -lm)
+  list(APPEND EXTERNAL_LIBS gptlf gptl)
 endif ()
 
 if (DEFINED ENV{PIO} AND (NOT DEFINED ENV{PIO_ROOT}))
@@ -398,27 +375,7 @@ if (HAS_RRTMGP AND HAS_NOAHMP)
   add_definitions(-DHAS_WRF)
 endif ()
 
-#find_package(OpenMP)
-#if(OpenMP_Fortran_FOUND)
-#  # 如果找到了OpenMP，为Fortran编译器设置OpenMP编译选项
-#  message("Found OMP!")
-#  set(OpenMP_Fortran_FLAGS "${OpenMP_Fortran_FLAGS} -qopenmp-link=static")
-#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
-#  message("${OpenMP_Fortran_LIB_NAMES}")
-#  message("${OpenMP_Fortran_LIBRARIES}")
-#  list(APPEND EXTERNAL_LIBS OpenMP::OpenMP_Fortran)
-#endif()
-
-if (DEFINED ENV{OMPDIR})
-  include_directories("$ENV{OMPDIR}/include")
-  link_directories("$ENV{OMPDIR}/lib")
-endif ()
-
-set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fiopenmp -qopenmp-link=static")
-
 target_link_libraries(gmcore fortran_container fortran_datetime fiona flogger ${EXTERNAL_LIBS})
-
-
 if (USE_CAM AND HAS_GPTL AND HAS_PIO AND HAS_LAPACK AND HAS_ESMF AND HAS_MCT)
   message(STATUS "GPTL, PIO, LAPACK, ESMF, and MCT are found. Building CAM physics and chemistry.")
   include_directories(src/physics/cam/include)
```

**src/drivers/gmcore_adv_driver.F90**
```fortran
diff --git a/gmcore/src/drivers/gmcore_adv_driver.F90 b/gmcore/src/drivers/gmcore_adv_driver.F90
index 9be9193..80f1a2d 100644
--- a/gmcore/src/drivers/gmcore_adv_driver.F90
+++ b/gmcore/src/drivers/gmcore_adv_driver.F90
@@ -20,7 +20,6 @@ program gmcore_adv_driver
   use deform_test_mod
   use moving_vortices_test_mod
   use dcmip12_test_mod
-  use perf_mod
 
   implicit none
 
@@ -49,7 +48,6 @@ program gmcore_adv_driver
   call gmcore_init_stage0()
 
   call gmcore_init_stage1(namelist_path)
-  call t_startf ( 'gmcore_inits' )
 
   select case (test_case)
   case ('solid_rotation')
@@ -94,34 +92,22 @@ program gmcore_adv_driver
   call history_setup_h0_adv(blocks)
   call output(old)
   call diagnose(old)
-
-  call t_stopf ( 'gmcore_inits' )
-
   if (proc%is_root()) call log_print_diag(curr_time%isoformat())
 
-  call t_startf ( 'run_not_gmcore' )
   time1 = MPI_WTIME()
   do while (.not. time_is_finished())
-    call t_startf ( 'loop_seg_1' )
     call set_uv(elapsed_seconds + dt_adv, new)
     call time_advance(dt_adv)
     call adv_accum_wind(old)
-    call t_stopf ( 'loop_seg_1' )
-
-    call t_startf ( 'adv_run' )
     call adv_run(old)
-    call t_stopf ( 'adv_run' )
     call diagnose(old)
     if (proc%is_root() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
     call output(old)
   end do
   time2 = MPI_WTIME()
   if (proc%is_root()) call log_notice('Total time cost ' // to_str(time2 - time1, 5) // ' seconds.')
-  call t_stopf ( 'run_not_gmcore' )
 
-  ! call t_startf ( 'gmcore_final' )
   call gmcore_final()
-  ! call t_stopf ( 'gmcore_final' )
 
 contains
```

**gmcore/src/drivers/gmcore_driver.F90**
```fortran
diff --git a/gmcore/src/drivers/gmcore_driver.F90 b/gmcore/src/drivers/gmcore_driver.F90
index bf144eb..1b156ef 100644
--- a/gmcore/src/drivers/gmcore_driver.F90
+++ b/gmcore/src/drivers/gmcore_driver.F90
@@ -28,7 +28,6 @@ program gmcore_driver
   use mars_cold_run_mod
   use tropical_cyclone_test_mod
   use prepare_mod
-  use perf_mod
 
   implicit none
 
@@ -62,7 +61,6 @@ program gmcore_driver
   end select
 
   call gmcore_init_stage1(namelist_path)
-  call t_startf ( 'gmcore_inits' )
 
   if (initial_file == 'N/A' .and. test_case == 'N/A' .and. .not. restart) then
     call prepare_topo()
@@ -115,14 +113,9 @@ program gmcore_driver
   end if
 
   call gmcore_init_stage3()
-  call t_stopf ( 'gmcore_inits' )
 
-  call t_startf ( 'gmcore_run' )
   call gmcore_run()
-  call t_stopf  ( 'gmcore_run' )
 
-  ! call t_startf ( 'gmcore_final' )
   call gmcore_final()
-  ! call t_stopf  ( 'gmcore_final' )
 
 end program gmcore_driver
```
**gmcore/src/drivers/gmcore_swm_driver.F90**
```fortran
diff --git a/gmcore/src/drivers/gmcore_swm_driver.F90 b/gmcore/src/drivers/gmcore_swm_driver.F90
index 436ec0e..ab753fa 100644
--- a/gmcore/src/drivers/gmcore_swm_driver.F90
+++ b/gmcore/src/drivers/gmcore_swm_driver.F90
@@ -1,7 +1,6 @@
 program gmcore_swm_driver
 
   use flogger
-  use string
   use const_mod
   use namelist_mod
   use block_mod
@@ -16,7 +15,6 @@ program gmcore_swm_driver
   use shallow_water_waves_test_mod
   use vortex_erosion_test_mod
   use splash_test_mod
-  use perf_mod
 
   implicit none
 
@@ -36,13 +34,10 @@ program gmcore_swm_driver
     call log_error('You should give a namelist file path!')
   end if
 
-
-
   call parse_namelist(namelist_path)
 
   call gmcore_init_stage0()
 
-  call t_startf ( 'gmcore_inits' )
   call gmcore_init_stage1(namelist_path)
 
   select case (test_case)
@@ -52,8 +47,6 @@ program gmcore_swm_driver
 
   call gmcore_init_stage2(namelist_path)
 
-  call t_stopf  ( 'gmcore_inits' )
-
   if (restart) then
     call restart_read()
   else
@@ -83,15 +76,8 @@ program gmcore_swm_driver
     end do
   end if
 
-    ! print *, "Size of blocks array:", size(blocks)
-
-  call t_startf ( 'gmcore_run' )
   call gmcore_run()
-  ! call exbdrift_ion_vels( lchnk, ncol, pbuf)
-  call t_stopf  ( 'gmcore_run' )
 
-  ! call t_startf ( 'gmcore_final' )
   call gmcore_final()
-  ! call t_startf ( 'gmcore_final' )
 
 end program gmcore_swm_driver
```

**gmcore/src/dynamics/adv/adv_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/adv/adv_mod.F90 b/gmcore/src/dynamics/adv/adv_mod.F90
index 374eeb3..0f8d890 100644
--- a/gmcore/src/dynamics/adv/adv_mod.F90
+++ b/gmcore/src/dynamics/adv/adv_mod.F90
@@ -26,7 +26,6 @@ module adv_mod
   use weno_mod
   use tvd_mod
   use physics_mod
-  use perf_mod
 
   implicit none
 
@@ -135,14 +134,12 @@ contains
     type(latlon_field3d_type), intent(inout) :: qmfy
     real(r8), intent(in), optional :: dt
 
-    call perf_start('adv_calc_tracer_hflx') 
     select case (batch%scheme)
     case ('upwind')
       call upwind_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)
     case ('ffsl')
       call ffsl_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)
     end select
-    call perf_stop('adv_calc_tracer_hflx') 
 
   end subroutine adv_calc_tracer_hflx
 
@@ -152,15 +149,13 @@ contains
     type(latlon_field3d_type), intent(in   ) :: q
     type(latlon_field3d_type), intent(inout) :: qmfz
     real(r8), intent(in), optional :: dt
-    
-    call perf_start('adv_calc_tracer_vflx') 
+
     select case (batch%scheme)
     case ('upwind')
       call upwind_calc_tracer_vflx(batch, q, qmfz, dt)
     case ('ffsl')
       call ffsl_calc_tracer_vflx(batch, q, qmfz, dt)
     end select
-    call perf_stop('adv_calc_tracer_vflx') 
 
   end subroutine adv_calc_tracer_vflx
```
**gmcore/src/dynamics/adv/ffsl_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/adv/ffsl_mod.F90 b/gmcore/src/dynamics/adv/ffsl_mod.F90
index 10c6ee7..8b6e95d 100644
--- a/gmcore/src/dynamics/adv/ffsl_mod.F90
+++ b/gmcore/src/dynamics/adv/ffsl_mod.F90
@@ -27,7 +27,7 @@ module ffsl_mod
   use adv_batch_mod
   use ppm_mod
   use limiter_mod
-  use perf_mod
+
   implicit none
 
   private
@@ -93,8 +93,6 @@ contains
     real(r8) pole(m%mesh%full_nlev)
     real(r8) dt_opt
 
-    call perf_start('ffsl_calc_mass_hflx')
-
     dt_opt = batch%dt; if (present(dt)) dt_opt = dt
 
     associate (mesh => m%mesh    , &
@@ -168,8 +166,6 @@ contains
     call hflx(batch, u, v, my, mx, mfx, mfy)
     end associate
 
-    call perf_stop('ffsl_calc_mass_hflx')
-
   end subroutine ffsl_calc_mass_hflx
 
   subroutine ffsl_calc_mass_vflx(batch, m, mfz, dt)
@@ -178,9 +174,9 @@ contains
     type(latlon_field3d_type), intent(in   ) :: m
     type(latlon_field3d_type), intent(inout) :: mfz
     real(r8), intent(in), optional :: dt
-    call perf_start('ffsl_calc_mass_vflx')
+
     call vflx(batch, batch%we, m, mfz)
-    call perf_stop('ffsl_calc_mass_vflx')
+
   end subroutine ffsl_calc_mass_vflx
 
   subroutine ffsl_calc_tracer_hflx(batch, q, qmfx, qmfy, dt)
@@ -196,8 +192,6 @@ contains
     real(r8) pole(q%mesh%half_nlev)
     real(r8) dt_opt
 
-    call perf_start('ffsl_calc_tracer_hflx')
-
     dt_opt = batch%dt; if (present(dt)) dt_opt = dt
 
     associate (mesh => q%mesh    , &
@@ -281,7 +275,7 @@ contains
     ! Run outer flux form operators.
     call hflx(batch, mfx, mfy, qy, qx, qmfx, qmfy)
     end associate
-    call perf_stop('ffsl_calc_tracer_hflx')
+
   end subroutine ffsl_calc_tracer_hflx
 
   subroutine ffsl_calc_tracer_vflx(batch, q, qmfz, dt)
@@ -290,11 +284,9 @@ contains
     type(latlon_field3d_type), intent(in   ) :: q
     type(latlon_field3d_type), intent(inout) :: qmfz
     real(r8), intent(in), optional :: dt
-    
-    call perf_start('ffsl_calc_tracer_vflx')
+
     call vflx(batch, batch%we, q, qmfz)
 
-    call perf_stop('ffsl_calc_tracer_vflx')
   end subroutine ffsl_calc_tracer_vflx
 
   subroutine hflx_van_leer(batch, u, v, mx, my, mfx, mfy)
@@ -310,7 +302,6 @@ contains
     integer ks, ke, i, j, k, iu, ju, ci
     real(r8) cf, dm
 
-    call perf_start('hflx_van_leer')
     associate (mesh => u%mesh    , &
                cflx => batch%cflx, & ! in
                cfly => batch%cfly)   ! in
@@ -351,8 +342,6 @@ contains
     end select
     end associate
 
-    call perf_stop('hflx_van_leer')
-
   end subroutine hflx_van_leer
 
   subroutine vflx_van_leer(batch, w, m, mfz)
@@ -365,9 +354,6 @@ contains
     integer i, j, k, ku, ci
     real(r8) cf, dm
 
-
-    call perf_start('vflx_van_leer')
-
     associate (mesh => m%mesh    , &
                cflz => batch%cflz)   ! in
     select case (batch%loc)
@@ -414,8 +400,6 @@ contains
     end select
     end associate
 
-    call perf_stop('vflx_van_leer')
-
   end subroutine vflx_van_leer
 
   subroutine hflx_ppm(batch, u, v, mx, my, mfx, mfy)
@@ -431,7 +415,6 @@ contains
     integer ks, ke, i, j, k, iu, ju, ci
     real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6
 
-    call perf_start('hflx_ppm')
     associate (mesh => u%mesh    , &
                cflx => batch%cflx, & ! in
                cfly => batch%cfly)   ! in
@@ -439,8 +422,6 @@ contains
     case ('cell', 'lev')
       ks = merge(mesh%full_kds, mesh%half_kds, batch%loc == 'cell')
       ke = merge(mesh%full_kde, mesh%half_kde, batch%loc == 'cell')
-     !$omp parallel 
-     !$omp do private(i, j, k, iu, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(2)
       do k = ks, ke
         ! Along x-axis
         do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
@@ -470,10 +451,6 @@ contains
             end if
           end do
         end do
-      end do 
-      !$omp end do
-      !$omp do private(i, j, k, ju, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(2)
-      do k = ks, ke
         ! Along y-axis
         do j = mesh%half_jds, mesh%half_jde
           do i = mesh%full_ids, mesh%full_ide
@@ -501,14 +478,10 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel 
     case ('vtx')
     end select
     end associate
 
-    call perf_stop('hflx_ppm')
-  
   end subroutine hflx_ppm
 
   subroutine vflx_ppm(batch, w, m, mfz)
@@ -521,13 +494,10 @@ contains
     integer i, j, k, ku, ci
     real(r8) cf, s1, s2, ds1, ds2, ds3, ml, dm, m6
 
-    call perf_start('vflx_ppm')
-    
     associate (mesh => m%mesh    , &
                cflz => batch%cflz)   ! in
     select case (batch%loc)
     case ('cell')
-      !$omp parallel do private(i, j, k, ku, ml, dm, m6, s1, s2, ds1, ds2, ds3, cf, ci) collapse(2)
       do k = mesh%half_kds + 1, mesh%half_kde - 1
         do j = mesh%full_jds, mesh%full_jde
           do i = mesh%full_ids, mesh%full_ide
@@ -559,9 +529,7 @@ contains
           end do
         end do
       end do
-      !$omp end parallel do
     case ('lev')
-      !$omp parallel do private(k, j, i, ci, cf, ku, ml, dm, m6, s1, s2, ds1, ds2, ds3) collapse(2)
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds, mesh%full_jde
           do i = mesh%full_ids, mesh%full_ide
@@ -591,12 +559,9 @@ contains
           end do
         end do
       end do
-      !$omp end parallel do
     end select
     end associate
 
-    call perf_stop('vflx_ppm')
-
   end subroutine vflx_ppm
 
 end module ffsl_mod
```
**gmcore/src/dynamics/damp/damp_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/damp/damp_mod.F90 b/gmcore/src/dynamics/damp/damp_mod.F90
index 212f7c7..9e871d0 100644
--- a/gmcore/src/dynamics/damp/damp_mod.F90
+++ b/gmcore/src/dynamics/damp/damp_mod.F90
@@ -16,7 +16,6 @@ module damp_mod
   use vor_damp_mod
   use smag_damp_mod
   use laplace_damp_mod
-  use perf_mod
 
   implicit none
 
@@ -56,7 +55,6 @@ contains
 
     ! This nudging of polar v helps to keep the flow neat around the poles.
     ! NOTE: DO NOT REMOVE IT!
-    call t_startf ( 'damp_run' )
     do j = block%mesh%half_jms, block%mesh%half_jme
       if (block%mesh%is_south_pole(j)) then
         dstate%v_lat%d(:,j,:) = 0.8_r8 * dstate%v_lat%d(:,j,:) + 0.2_r8 * dstate%v_lat%d(:,j+1,:)
@@ -80,7 +78,6 @@ contains
         call smag_damp_run(block, dstate, dt / smag_damp_cycles)
       end do
     end if
-    call t_stopf ( 'damp_run' )
 
   end subroutine damp_run
```
**gmcore/src/dynamics/damp/div_damp_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/damp/div_damp_mod.F90 b/gmcore/src/dynamics/damp/div_damp_mod.F90
index 3af3247..b446e77 100644
--- a/gmcore/src/dynamics/damp/div_damp_mod.F90
+++ b/gmcore/src/dynamics/damp/div_damp_mod.F90
@@ -90,7 +90,7 @@ contains
   end subroutine div_damp_final
 
   subroutine div_damp_run(block, dstate, dt)
-    ! Add omp here cause errors
+
     type(block_type), intent(inout) :: block
     type(dstate_type), intent(inout) :: dstate
     real(8), intent(in) :: dt
```

**gmcore/src/dynamics/filter_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/filter_mod.F90 b/gmcore/src/dynamics/filter_mod.F90
index 24b9a2d..3b256d2 100644
--- a/gmcore/src/dynamics/filter_mod.F90
+++ b/gmcore/src/dynamics/filter_mod.F90
@@ -112,8 +112,6 @@ contains
       ngrid => filter%ngrid_lon
     end select
 
-    !$omp parallel
-    !$omp do collapse(2) private(i, j, k, n, hn, tmp)
     do k = ks, ke
       do j = js, je
         if (ngrid(j) > 1) then
@@ -132,8 +130,6 @@ contains
         end if
       end do
     end do
-    !$omp end do
-    !$omp end parallel
 
     call perf_stop('filter_run_3d')
```

**gmcore/src/dynamics/gmcore_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/gmcore_mod.F90 b/gmcore/src/dynamics/gmcore_mod.F90
index 2c969fc..f92a95f 100644
--- a/gmcore/src/dynamics/gmcore_mod.F90
+++ b/gmcore/src/dynamics/gmcore_mod.F90
@@ -137,7 +137,6 @@ contains
 
     call global_mesh%init_global(nlon, nlat, nlev, lon_hw=lon_hw, lat_hw=lat_hw)
     call process_create_blocks()
-    ! print *, "Size of blocks array:", size(blocks)
     associate (mesh => blocks(1)%mesh)
     min_lon = mesh%full_lon_deg(mesh%full_ims)
     max_lon = mesh%full_lon_deg(mesh%full_ime)
@@ -199,8 +198,6 @@ contains
 
     if (proc%is_root()) call print_namelist()
 
-    ! print *, "Size of blocks array:", size(blocks)
-
     do iblk = 1, size(blocks)
       blocks(iblk)%mesh%full_lev  = global_mesh%full_lev
       blocks(iblk)%mesh%half_lev  = global_mesh%half_lev
@@ -227,8 +224,6 @@ contains
   subroutine gmcore_run()
 
     integer i, j, m, iblk, itime
-    ! perf_start()
-    call t_startf ( 'aaagmcore_run_inits' )
 
     do iblk = 1, size(blocks)
       associate (block => blocks(iblk)     , &
@@ -257,20 +252,12 @@ contains
     if (proc%is_root()) call log_print_diag(curr_time%isoformat())
     call output(old)
 
-    call t_stopf ( 'aaagmcore_run_inits' )
-
-    call t_startf ( 'aaagmcore_main_loop' )
-
     model_main_loop: do while (.not. time_is_finished())
       ! ------------------------------------------------------------------------
       !                              Dynamical Core
-      call t_startf ( 'aaagmcore_dynamic_core' )
       do iblk = 1, size(blocks)
-        call t_startf ( 'time_integrator' )
         call time_integrator(operators, blocks(iblk), old, new, dt_dyn)
-        call t_stopf ( 'time_integrator' )
         call damp_run(blocks(iblk), blocks(iblk)%dstate(new), dt_dyn)
-        ! call t_stopf ( 'aaagmcore_damp_run' )
         if (pdc_type == 1) call physics_update_dynamics(blocks(iblk), new, dt_dyn)
         call blocks(iblk)%dstate(new)%c2a()
       end do
@@ -278,15 +265,11 @@ contains
       ! Advance to n+1 time level.
       ! NOTE: Time indices are swapped, e.g. new <=> old.
       call time_advance(dt_dyn)
-      call t_stopf ( 'aaagmcore_dynamic_core' )
       ! ------------------------------------------------------------------------
       !                            Tracer Advection
-      call t_startf ( 'aaagmcore_adv_run' )
       call adv_run(old)
-      call t_stopf ( 'aaagmcore_adv_run' )
       ! ------------------------------------------------------------------------
       !                                Physics
-      call t_startf ( 'aaagmcore_physic_run' )
       call test_forcing_run(dt_dyn, old)
       if (baroclinic) then
         do iblk = 1, size(blocks)
@@ -294,7 +277,6 @@ contains
           if (pdc_type == 3) call physics_update(blocks(iblk), old, dt_phys)
         end do
       end if
-      call t_stopf ( 'aaagmcore_physic_run' )
       ! ------------------------------------------------------------------------
       call diagnose(blocks, old)
       if (proc%is_root() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
@@ -302,11 +284,10 @@ contains
       call output(old)
     end do model_main_loop
 
-    call t_stopf ( 'aaagmcore_main_loop' )
-
   end subroutine gmcore_run
 
   subroutine gmcore_final()
+
     call log_final()
     call time_final()
     call interp_final()
@@ -495,8 +476,6 @@ contains
 
     integer i, j, k
 
-    call t_startf ('space_operators')
-
     call dtend1%reset_flags()
 
     associate (mesh => block%mesh)
@@ -618,9 +597,6 @@ contains
     end select
     end associate
 
-
-    call t_stopf ('space_operators')
-
   end subroutine space_operators
 
 end module gmcore_mod
```

**gmcore/src/dynamics/interp_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/interp_mod.F90 b/gmcore/src/dynamics/interp_mod.F90
index b07167c..0d50bed 100644
--- a/gmcore/src/dynamics/interp_mod.F90
+++ b/gmcore/src/dynamics/interp_mod.F90
@@ -13,7 +13,6 @@ module interp_mod
   use namelist_mod
   use latlon_mesh_mod, mesh_type => latlon_mesh_type
   use latlon_field_types_mod
-  use perf_mod
 
   implicit none
 
@@ -87,8 +86,6 @@ contains
     select case (trim(x%loc) // '>' // trim(y%loc))
     ! --------------------------------------------------------------------------
     case ('cell>lon')
-    !$omp parallel 
-    !$omp do private(i, j, k) collapse(2) 
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
           do i = x%mesh%half_ids, x%mesh%half_ide
@@ -98,12 +95,8 @@ contains
           end do
         end do
       end do
-    !$omp end do
-    !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('cell>lat')
-    !$omp parallel 
-    !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%half_jds, x%mesh%half_jde
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -113,8 +106,6 @@ contains
           end do
         end do
       end do
-    !$omp end do
-    !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('cell>lev')
       if (x%mesh%full_nlev == 1) return
@@ -127,8 +118,6 @@ contains
       ! ===o=== k
       !
       ! -------
-    !$omp parallel 
-    !$omp do collapse(1) private(i, j, k, a, b, x1, x2)
       do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
         a = x%mesh%full_dlev(k-1) / (2 * x%mesh%half_dlev(k))
         b = x%mesh%full_dlev(k  ) / (2 * x%mesh%half_dlev(k))
@@ -138,8 +127,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
       k = x%mesh%half_kds
       x1 = x%mesh%full_lev(k  ) - x%mesh%half_lev(k)
       x2 = x%mesh%full_lev(k+1) - x%mesh%half_lev(k)
@@ -160,11 +147,8 @@ contains
           y%d(i,j,k) = a * x%d(i,j,k-1) + b * x%d(i,j,k-2)
         end do
       end do
-
     ! --------------------------------------------------------------------------
     case ('cell>vtx')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%half_jds, x%mesh%half_jde
           do i = x%mesh%half_ids, x%mesh%half_ide
@@ -175,12 +159,8 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('lon>cell')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -190,12 +170,8 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('lat>cell')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -205,8 +181,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('lev>cell')
       ! =======
@@ -218,8 +192,6 @@ contains
       ! ---o--- k+1
       !
       ! =======
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%full_jds, x%mesh%full_jde
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -227,8 +199,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('lon>lev_lon')
       ! -------
@@ -240,8 +210,6 @@ contains
       ! ===o=== k
       !
       ! ----o--
-      !$omp parallel 
-      !$omp do collapse(1) private(i, j, k, a, b)
       do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
         a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
         b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
@@ -251,9 +219,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
-
       k = x%mesh%half_kds
       ! ---?--- 1
       !
@@ -311,8 +276,6 @@ contains
       ! ===o=== k
       !
       ! -------
-      !$omp parallel 
-      !$omp do collapse(1) private(i, j, k, a, b)
       do k = x%mesh%half_kds + 1, x%mesh%half_kde - 1
         a = x%mesh%full_dlev(k-1) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
         b = x%mesh%full_dlev(k  ) / (x%mesh%full_dlev(k-1) + x%mesh%full_dlev(k))
@@ -322,8 +285,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
       k = x%mesh%half_kds
       ! ---?--- 1
       !
@@ -372,8 +333,6 @@ contains
       end do
     ! --------------------------------------------------------------------------
     case ('lev>lev_lon')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%half_kds, x%mesh%half_kde
         do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
           do i = x%mesh%half_ids, x%mesh%half_ide
@@ -383,12 +342,8 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('lev>lev_lat')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%half_kds, x%mesh%half_kde
         do j = x%mesh%half_jds, x%mesh%half_jde
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -398,8 +353,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     end select
 
   end subroutine interp_run_3d
@@ -411,15 +364,10 @@ contains
 
     integer i, j, k
 
-    ! real(r8) tempx(x%mesh%full_jds:x%mesh%full_jde+1, x%mesh%full_ids:x%mesh%full_ide, x%mesh%full_kds:x%mesh%full_kde)
-
-    call t_startf ('average_run')
 
     select case (trim(x%loc) // '>' // trim(y%loc))
     ! --------------------------------------------------------------------------
     case ('cell>lon')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%full_jds_no_pole, x%mesh%full_jde_no_pole
           do i = x%mesh%half_ids, x%mesh%half_ide
@@ -427,12 +375,8 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     ! --------------------------------------------------------------------------
     case ('cell>lat')
-      !$omp parallel 
-      !$omp do collapse(2) private(i, j, k)
       do k = x%mesh%full_kds, x%mesh%full_kde
         do j = x%mesh%half_jds, x%mesh%half_jde
           do i = x%mesh%full_ids, x%mesh%full_ide
@@ -440,8 +384,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     end select
 
   end subroutine average_run_3d
```

**gmcore/src/dynamics/nh_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/nh_mod.F90 b/gmcore/src/dynamics/nh_mod.F90
index e4aa230..aff6527 100644
--- a/gmcore/src/dynamics/nh_mod.F90
+++ b/gmcore/src/dynamics/nh_mod.F90
@@ -21,7 +21,6 @@ module nh_mod
   use tracer_mod
   use operators_mod
   use filter_mod
-  use perf_mod
 
   implicit none
 
@@ -39,7 +38,6 @@ contains
     type(dstate_type), intent(inout) :: new_dstate
     real(r8), intent(in) :: dt
 
-    call t_startf ('nh_solve')
     call interp_wind(block, star_dstate)
     call calc_adv_lev(block, star_dstate%w_lev , block%aux%adv_w_lev , star_dstate%dmg_lev, dt)
     call calc_adv_lev(block, star_dstate%gz_lev, block%aux%adv_gz_lev, star_dstate%dmg_lev, dt)
@@ -48,7 +46,6 @@ contains
     call calc_p(block, new_dstate)
     call interp_run(new_dstate%gz_lev, new_dstate%gz)
     call fill_halo(new_dstate%gz, west_halo=.false., south_halo=.false.)
-    call t_stopf ('nh_solve')
 
   end subroutine nh_solve
```

**gmcore/src/dynamics/operators_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/operators_mod.F90 b/gmcore/src/dynamics/operators_mod.F90
index 0f2403e..c50b802 100644
--- a/gmcore/src/dynamics/operators_mod.F90
+++ b/gmcore/src/dynamics/operators_mod.F90
@@ -14,7 +14,6 @@ module operators_mod
   use adv_mod
   use interp_mod
   use filter_mod
-  use omp_lib
 
   implicit none
 
@@ -147,9 +146,6 @@ contains
                mgs     => dstate%mgs    , & ! in
                mg_lev  => dstate%mg_lev , & ! out
                mg      => dstate%mg     )   ! out
-    
-    !$omp parallel 
-    !$omp do private(i, j, k) collapse(2)
     do k = mesh%half_kds, mesh%half_kde
       do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -157,8 +153,6 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp do private(i, j, k) collapse(2)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -166,8 +160,6 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp end parallel 
     end associate
 
     call perf_stop('calc_mg')
@@ -175,7 +167,7 @@ contains
   end subroutine calc_mg
 
   subroutine calc_ph(block, dstate)
-  ! Attention! Dependency!
+
     type(block_type), intent(inout) :: block
     type(dstate_type), intent(inout) :: dstate
 
@@ -203,8 +195,6 @@ contains
         end do
       end do
     end do
-    !$omp parallel 
-    !$omp do collapse(2) private(i, j, k)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -212,8 +202,6 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp end parallel 
     ! NOTE: Move this to other place?
     if (hydrostatic) ps%d = phs%d
     end associate
@@ -304,7 +292,6 @@ contains
                t    => dstate%t           , & ! out
                tv   => dstate%tv          )   ! out
     if (idx_qv > 0) then
-      !$omp parallel do collapse(2) private(k, j, i)
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
           do i = mesh%full_ids, mesh%full_ide + 1
@@ -313,9 +300,7 @@ contains
           end do
         end do
       end do
-      !$omp end parallel do
     else
-      !$omp parallel do collapse(2) private(k, j, i)
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
           do i = mesh%full_ids, mesh%full_ide + 1
@@ -324,7 +309,6 @@ contains
           end do
         end do
       end do
-      !$omp end parallel do
     end if
     end associate
 
@@ -366,9 +350,6 @@ contains
     real(r8), intent(in) :: dt
 
     integer i, j, k
-    integer result
-    real(r8) sum_dmf(block%mesh%full_ids:block%mesh%full_ide, &
-                     block%mesh%full_jds:block%mesh%full_jde)
 
     call perf_start('calc_we_lev')
 
@@ -378,18 +359,10 @@ contains
                we_lev     => dstate%we_lev       , & ! out
                we_lev_lon => block%aux%we_lev_lon, & ! out
                we_lev_lat => block%aux%we_lev_lat)   ! out
-
-
-    do j = mesh%full_jds, mesh%full_jde
-      do i = mesh%full_ids, mesh%full_ide
-        do k = mesh%half_kds + 1, mesh%half_kde - 1
-          if (k .eq. mesh%half_kds + 1) then 
-            sum_dmf(i,j) = sum(dmf%d(i,j,1:k-1))
-          else 
-            sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k-1)
-          end if
-          ! we_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum(dmf%d(i,j,1:k-1))
-          we_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum_dmf(i,j)
+    do k = mesh%half_kds + 1, mesh%half_kde - 1
+      do j = mesh%full_jds, mesh%full_jde
+        do i = mesh%full_ids, mesh%full_ide
+          we_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgs%d(i,j)) - sum(dmf%d(i,j,1:k-1))
         end do
       end do
     end do
@@ -409,7 +382,7 @@ contains
     type(dstate_type), intent(inout) :: dstate
 
     integer i, j, k
-    real(r8) ke_vtx_1, ke_vtx_2, ke_vtx_3, ke_vtx_4
+    real(r8) ke_vtx(4)
     real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
     real(r8) pole(block%mesh%full_nlev)
 
@@ -419,7 +392,6 @@ contains
                u    => dstate%u_lon, & ! in
                v    => dstate%v_lat, & ! in
                ke   => block%aux%ke)   ! out
-    !$omp parallel do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -431,7 +403,6 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
 
     if (ke_scheme == 2) then
       !
@@ -456,90 +427,72 @@ contains
       !     |________u________|________u________|
       !           i-1,j-1             i,j-1
       !
-      !$omp parallel do collapse(2) private(i, j, k, ke_vtx_1, ke_vtx_2, ke_vtx_3, ke_vtx_4)
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
           do i = mesh%full_ids, mesh%full_ide + 1
-            ke_vtx_1 = (                                    &
+            ke_vtx(1) = (                                    &
               mesh%area_lat_east (j  ) * v%d(i-1,j  ,k)**2 + &
               mesh%area_lat_west (j  ) * v%d(i  ,j  ,k)**2 + &
               mesh%area_lon_north(j  ) * u%d(i-1,j  ,k)**2 + &
               mesh%area_lon_south(j+1) * u%d(i-1,j+1,k)**2   &
             ) / mesh%area_vtx(j)
-            ke_vtx_2 = (                                    &
+            ke_vtx(2) = (                                    &
               mesh%area_lat_east (j-1) * v%d(i-1,j-1,k)**2 + &
               mesh%area_lat_west (j-1) * v%d(i  ,j-1,k)**2 + &
               mesh%area_lon_north(j-1) * u%d(i-1,j-1,k)**2 + &
               mesh%area_lon_south(j  ) * u%d(i-1,j  ,k)**2   &
             ) / mesh%area_vtx(j-1)
-            ke_vtx_3 = (                                    &
+            ke_vtx(3) = (                                    &
               mesh%area_lat_east (j-1) * v%d(i  ,j-1,k)**2 + &
               mesh%area_lat_west (j-1) * v%d(i+1,j-1,k)**2 + &
               mesh%area_lon_north(j-1) * u%d(i  ,j-1,k)**2 + &
               mesh%area_lon_south(j  ) * u%d(i  ,j  ,k)**2   &
             ) / mesh%area_vtx(j-1)
-            ke_vtx_4 = (                                    &
+            ke_vtx(4) = (                                    &
               mesh%area_lat_east (j  ) * v%d(i  ,j  ,k)**2 + &
               mesh%area_lat_west (j  ) * v%d(i+1,j  ,k)**2 + &
               mesh%area_lon_north(j  ) * u%d(i  ,j  ,k)**2 + &
               mesh%area_lon_south(j+1) * u%d(i  ,j+1,k)**2   &
             ) / mesh%area_vtx(j)
             ke%d(i,j,k) = (1.0_r8 - ke_cell_wgt) * (             &
-              (ke_vtx_1 + ke_vtx_4) * mesh%area_subcell(2,j) + &
-              (ke_vtx_2 + ke_vtx_3) * mesh%area_subcell(1,j)   &
+              (ke_vtx(1) + ke_vtx(4)) * mesh%area_subcell(2,j) + &
+              (ke_vtx(2) + ke_vtx(3)) * mesh%area_subcell(1,j)   &
             ) / mesh%area_cell(j) + ke_cell_wgt * ke%d(i,j,k)
           end do
         end do
       end do
-      !$omp end parallel do
     end if
 
     ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
     if (mesh%has_south_pole()) then
       j = mesh%full_jds
-      !!$omp parallel 
-      !!$omp do collapse(1) private(k, i)
       do k = mesh%full_kds, mesh%full_kde
         do i = mesh%full_ids, mesh%full_ide
           work(i,k) = 0.5_r8 * (v%d(i,j,k)**2 + block%aux%u_lat%d(i,j,k)**2)
         end do
       end do
-      !!$omp end do
-      !!$omp master
       call zonal_sum(proc%zonal_circle, work, pole)
       pole = pole / global_mesh%full_nlon
-      !!$omp end master
-      !!$omp do collapse(1) private(k, i)
       do k = mesh%full_kds, mesh%full_kde
         do i = mesh%full_ids, mesh%full_ide
           ke%d(i,j,k) = pole(k)
         end do
       end do
-      !!$omp end do
-      !!$omp end parallel
     end if
     if (mesh%has_north_pole()) then
       j = mesh%full_jde
-      !!$omp parallel 
-      !!$omp do collapse(1) private(k, i)
       do k = mesh%full_kds, mesh%full_kde
         do i = mesh%full_ids, mesh%full_ide
           work(i,k) = 0.5_r8 * (v%d(i,j-1,k)**2 + block%aux%u_lat%d(i,j-1,k)**2)
         end do
       end do
-      !!$omp end do
-      !!$omp master
       call zonal_sum(proc%zonal_circle, work, pole)
       pole = pole / global_mesh%full_nlon
-      !!$omp end master
-      !!$omp do collapse(1) private(k ,i)
       do k = mesh%full_kds, mesh%full_kde
         do i = mesh%full_ids, mesh%full_ide
           ke%d(i,j,k) = pole(k)
         end do
       end do
-      !!$omp end do
-      !!$omp end parallel
     end if
     end associate
 
@@ -639,18 +592,14 @@ contains
                ph_lev => dstate%ph_lev   , & ! in
                gz_lev => dstate%gz_lev   , & ! out
                gz     => dstate%gz       )   ! out
-    
     do k = mesh%half_kde - 1, mesh%half_kds, -1
-      !$omp parallel do collapse(1) private(i, j) 
       do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
           gz_lev%d(i,j,k) = gz_lev%d(i,j,k+1) + rd * tv%d(i,j,k) * log(ph_lev%d(i,j,k+1) / ph_lev%d(i,j,k))
         end do
       end do
-      !$omp end parallel do
     end do
     ! For output
-    !$omp parallel do collapse(2)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -658,67 +607,18 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
     end associate
 
     call perf_stop('calc_gz_lev')
 
   end subroutine calc_gz_lev
 
-  subroutine calc_gz_lev_rhod(block, dstate)
-
-    type(block_type), intent(in) :: block
-    type(dstate_type), intent(inout) :: dstate
-
-    integer i, j, k
-
-    call perf_start('calc_gz_lev_rhod')
-
-    associate (mesh   => block%mesh      , & ! in
-               gzs    => block%static%gzs, & ! in
-               tv     => dstate%tv       , & ! in
-               ph_lev => dstate%ph_lev   , & ! in
-               dmg    => dstate%dmg      , & ! in     
-               gz_lev => dstate%gz_lev   , & ! out
-               gz     => dstate%gz       , & ! out
-               rhod   => dstate%rhod  )   ! out
-    do k = mesh%half_kde - 1, mesh%half_kds, -1
-      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
-        do i = mesh%full_ids, mesh%full_ide + 1
-          gz_lev%d(i,j,k) = gz_lev%d(i,j,k+1) + rd * tv%d(i,j,k) * log(ph_lev%d(i,j,k+1) / ph_lev%d(i,j,k))
-        end do
-      end do
-    end do
-    ! For output
-    do k = mesh%full_kds, mesh%full_kde
-      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
-        do i = mesh%full_ids, mesh%full_ide + 1
-          gz%d(i,j,k) = 0.5_r8 * (gz_lev%d(i,j,k) + gz_lev%d(i,j,k+1))
-          rhod%d(i,j,k) = dmg%d(i,j,k) / (gz_lev%d(i,j,k) - gz_lev%d(i,j,k+1))
-        end do
-      end do
-    end do
-
-    ! do k = mesh%full_kds, mesh%full_kde
-    !   do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
-    !     do i = mesh%full_ids, mesh%full_ide + 1
-    !       rhod%d(i,j,k) = dmg%d(i,j,k) / (gz_lev%d(i,j,k) - gz_lev%d(i,j,k+1))
-    !     end do
-    !   end do
-    ! end do
-    end associate
-
-    call perf_stop('calc_gz_lev_rhod')
-
-  end subroutine calc_gz_lev_rhod
-
   subroutine calc_dmg(block, dstate)
 
     type(block_type), intent(inout) :: block
     type(dstate_type), intent(inout) :: dstate
 
     integer i, j, k, l
-    integer send_south_req, recv_south_req, send_north_req, recv_north_req
 
     call perf_start('calc_dmg')
 
@@ -733,8 +633,6 @@ contains
                dmg_lev => dstate%dmg_lev   , & ! out
                dmg_vtx => block%aux%dmg_vtx)   ! out
     if (baroclinic .or. advection) then
-      !$omp parallel 
-      !$omp do collapse(2)
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds, mesh%full_jde
           do i = mesh%full_ids, mesh%full_ide
@@ -751,8 +649,7 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp do collapse(2)
+
       do k = mesh%half_kds + 1, mesh%half_kde - 1
         do j = mesh%full_jds, mesh%full_jde
           do i = mesh%full_ids, mesh%full_ide
@@ -760,8 +657,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel 
       ! Top boundary
       k = mesh%half_kds
       do j = mesh%full_jds, mesh%full_jde
@@ -776,7 +671,7 @@ contains
           dmg_lev%d(i,j,k) = mg_lev%d(i,j,k) - mg%d(i,j,k-1)
         end do
       end do
-      ! call fill_halo(dmg_lev)
+      call fill_halo(dmg_lev)
     else
       do j = mesh%full_jds, mesh%full_jde
         do i = mesh%full_ids, mesh%full_ide
@@ -804,7 +699,6 @@ contains
     real(r8), intent(in) :: dt
 
     integer i, j, k
-    integer send_south_req, recv_south_req, send_north_req, recv_north_req
 
     call perf_start('calc_mf')
 
@@ -820,8 +714,6 @@ contains
                mfy_lat => block%aux%mfy_lat, & ! out
                mfy_lon => block%aux%mfy_lon, & ! out
                mfx_lat => block%aux%mfx_lat)   ! out
-    !$omp parallel 
-    !$omp do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
         do i = mesh%half_ids - 1, mesh%half_ide
@@ -829,8 +721,6 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
         do i = mesh%full_ids, mesh%full_ide + 1
@@ -838,8 +728,7 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp do collapse(2) private(k, j, i)
+
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%half_jds, mesh%half_jde
         do i = mesh%full_ids, mesh%full_ide
@@ -849,13 +738,8 @@ contains
         end do
       end do
     end do
-    !$omp end do
-
-    !$omp master
-    ! call fill_halo(u_lat)
-    !$omp end master
+    call fill_halo(u_lat)
 
-    !$omp do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
         do i = mesh%half_ids, mesh%half_ide
@@ -865,9 +749,7 @@ contains
         end do
       end do
     end do
-    !$omp end do
-    !$omp end parallel
-    ! call fill_halo(v_lon)
+    call fill_halo(v_lon)
     end associate
 
     call perf_stop('calc_mf')
@@ -890,7 +772,6 @@ contains
                v_lat => dstate%v_lat   , & ! in
                u_lat => block%aux%u_lat, & ! in
                vor   => block%aux%vor  )   ! out
-    !$omp parallel do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%half_jds, mesh%half_jde
         do i = mesh%half_ids, mesh%half_ide
@@ -901,57 +782,37 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
     if (pv_pole_stokes) then
       ! Special treatment of vorticity around Poles
       if (mesh%has_south_pole()) then
         j = mesh%half_jds
-        !
-        
-        !
-        !!$omp parallel 
-        !!$omp do collapse(2) private(k, i)
         do k = mesh%full_kds, mesh%full_kde
           do i = mesh%half_ids, mesh%half_ide
             work(i,k) = -u_lat%d(i,j,k) * mesh%le_lat(j)
           end do
         end do
-        !!$omp end do
-        !!$omp master
         call zonal_sum(proc%zonal_circle, work, pole)
         pole = pole / global_mesh%full_nlon / mesh%area_cell(j)
-        !!$omp end master
-        !!$omp do collapse(2) private(k, i)
         do k = mesh%full_kds, mesh%full_kde
           do i = mesh%half_ids, mesh%half_ide
             vor%d(i,j,k) = pole(k)
           end do
         end do
-        !!$omp end do
-        !!$omp end parallel
       end if
       if (mesh%has_north_pole()) then
         j = mesh%half_jde
-        !!$omp parallel
-        !!$omp do collapse(2) private(i, k)
         do k = mesh%full_kds, mesh%full_kde
           do i = mesh%half_ids, mesh%half_ide
             work(i,k) = u_lat%d(i,j,k) * mesh%le_lat(j)
           end do
         end do
-        !!$omp end do
-        !!$omp master
         call zonal_sum(proc%zonal_circle, work, pole)
         pole = pole / global_mesh%full_nlon / mesh%area_cell(j+1)
-        !!$omp end master
-        !!$omp do collapse(1) private(k, i)
         do k = mesh%full_kds, mesh%full_kde
           do i = mesh%half_ids, mesh%half_ide
             vor%d(i,j,k) = pole(k)
           end do
         end do
-        !!$omp end do
-        !!$omp end parallel
       end if
     end if
     end associate
@@ -974,7 +835,6 @@ contains
                vor     => block%aux%vor    , & ! in
                pv      => block%aux%pv     )   ! out
     call calc_vor(block, dstate)
-    !$omp parallel do collapse(2) private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%half_jds, mesh%half_jde
         do i = mesh%half_ids, mesh%half_ide
@@ -982,7 +842,6 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
     call fill_halo(pv)
     end associate
 
@@ -1050,8 +909,6 @@ contains
                pv_lat => block%aux%pv_lat)   ! out
     select case (upwind_order_pv)
     case (1)
-      !$omp parallel
-      !$omp do collapse(2) private(b, k, j, i) 
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
           do i = mesh%half_ids, mesh%half_ide
@@ -1060,10 +917,6 @@ contains
                               (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
           end do
         end do
-      end do
-      !$omp end do
-      !$omp do collapse(2) private(b, k, j, i) 
-      do k = mesh%full_kds, mesh%full_kde
         do j = mesh%half_jds, mesh%half_jde
           do i = mesh%full_ids, mesh%full_ide
             b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
@@ -1072,11 +925,7 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     case (3)
-      !$omp parallel
-      !$omp do collapse(2) private(b, k, j, i) 
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
           do i = mesh%half_ids, mesh%half_ide
@@ -1085,10 +934,6 @@ contains
                               (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
           end do
         end do
-      end do
-      !$omp end do
-      !$omp do collapse(2) private(b, k, j, i) 
-      do k = mesh%full_kds, mesh%full_kde
         do j = mesh%half_jds, mesh%half_jde
           do i = mesh%full_ids, mesh%full_ide
             b  = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
@@ -1097,11 +942,7 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     case (5)
-      !$omp parallel
-      !$omp do collapse(2) private(b, k, j, i) 
       do k = mesh%full_kds, mesh%full_kde
         do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
           do i = mesh%half_ids, mesh%half_ide
@@ -1110,11 +951,6 @@ contains
                               (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
           end do
         end do
-      end do
-      !$omp end do
-
-      !$omp do collapse(2) private(b, k, j, i)
-      do k = mesh%full_kds, mesh%full_kde
         do j = mesh%half_jds, mesh%half_jde
           do i = mesh%full_ids, mesh%full_ide
             b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
@@ -1123,8 +959,6 @@ contains
           end do
         end do
       end do
-      !$omp end do
-      !$omp end parallel
     end select
     call fill_halo(pv_lon, east_halo=.false., south_halo=.false.)
     call fill_halo(pv_lat, west_halo=.false., north_halo=.false.)
@@ -1357,7 +1191,6 @@ contains
     call adv_calc_tracer_hflx(block%adv_batch_pt, pt, ptf_lon, ptf_lat, dt)
     call fill_halo(ptf_lon, south_halo=.false., north_halo=.false., east_halo=.false.)
     call fill_halo(ptf_lat, north_halo=.false.,  west_halo=.false., east_halo=.false.)
-    !$omp parallel do private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
         do i = mesh%full_ids, mesh%full_ide
@@ -1370,7 +1203,6 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
     if (mesh%has_south_pole()) then
       j = mesh%full_jds
       do k = mesh%full_kds, mesh%full_kde
@@ -1404,7 +1236,6 @@ contains
     ! --------------------------------- FFSL -----------------------------------
     call adv_fill_vhalo(pt)
     call adv_calc_tracer_vflx(block%adv_batch_pt, pt, ptf_lev, dt)
-    !$omp parallel do private(k, j, i)
     do k = mesh%full_kds, mesh%full_kde
       do j = mesh%full_jds, mesh%full_jde
         do i = mesh%full_ids, mesh%full_ide
@@ -1412,7 +1243,6 @@ contains
         end do
       end do
     end do
-    !$omp end parallel do
     end associate
 
     call perf_stop('calc_grad_ptf')
```

**gmcore/src/dynamics/time_schemes_mod.F90**
```fortran
diff --git a/gmcore/src/dynamics/time_schemes_mod.F90 b/gmcore/src/dynamics/time_schemes_mod.F90
index 1758b40..26653c1 100644
--- a/gmcore/src/dynamics/time_schemes_mod.F90
+++ b/gmcore/src/dynamics/time_schemes_mod.F90
@@ -27,7 +27,6 @@ module time_schemes_mod
   use latlon_parallel_mod
   use process_mod, only: proc
   use filter_mod
-  use perf_mod
 
   implicit none
 
@@ -157,9 +156,6 @@ contains
 
     integer i, j, k
 
-    call t_startf ('update_state')
-
-
     associate (mesh       => block%mesh, &
                dmgsdt     => dtend%dmgs, &
                dgzdt      => dtend%dgz , &
@@ -252,8 +248,6 @@ contains
     end if
     end associate
 
-    call t_stopf ('update_state')
-
   end subroutine update_state
 
   subroutine rk2(space_operators, block, old, new, dt)
```

**gmcore/src/dynamics/time_schemes_mod.F90**
```fortran
diff --git a/gmcore/src/meshes/latlon/latlon_decomp_mod.F90 b/gmcore/src/meshes/latlon/latlon_decomp_mod.F90
index 9a75731..7ab5583 100644
--- a/gmcore/src/meshes/latlon/latlon_decomp_mod.F90
+++ b/gmcore/src/meshes/latlon/latlon_decomp_mod.F90
@@ -168,9 +168,6 @@ contains
     end if
 
     ! Set grid_proc_idmap for later use.
-    ! write(*, *) "nlon:", nlon
-    ! write(*, *) "nlat:", nlat
-
     allocate(proc%grid_proc_idmap(nlon,nlat))
     allocate(proc%global_grid_id (nlon,nlat))
     allocate(proc%local_grid_id  (nlon,nlat))
```

**gmcore/src/dynamics/time_schemes_mod.F90**
```fortran
diff --git a/gmcore/src/meshes/latlon/latlon_mesh_mod.F90 b/gmcore/src/meshes/latlon/latlon_mesh_mod.F90
index ddc9f6d..2c9be7f 100644
--- a/gmcore/src/meshes/latlon/latlon_mesh_mod.F90
+++ b/gmcore/src/meshes/latlon/latlon_mesh_mod.F90
@@ -150,11 +150,6 @@ contains
 
     call this%clear(keep_lev)
 
-
-    ! write(*, *) "Nlon:", nlon
-    ! write(*, *) "Nlat:", nlat
-    ! write(*, *) "Nlev:", nlev
-
     this%full_nlon = nlon
     this%half_nlon = nlon
     this%full_ids  = 1
```
**gmcore/src/meshes/latlon/latlon_parallel_mod.F90**
```fortran
diff --git a/gmcore/src/meshes/latlon/latlon_parallel_mod.F90 b/gmcore/src/meshes/latlon/latlon_parallel_mod.F90
index 01ab010..f18166f 100644
--- a/gmcore/src/meshes/latlon/latlon_parallel_mod.F90
+++ b/gmcore/src/meshes/latlon/latlon_parallel_mod.F90
@@ -34,8 +34,6 @@ module latlon_parallel_mod
   interface fill_halo
     module procedure fill_halo_2d
     module procedure fill_halo_3d
-    module procedure fill_halo_start_3d
-    module procedure fill_halo_stop_3d
     module procedure fill_halo_4d
   end interface fill_halo
 
@@ -180,38 +178,11 @@ contains
     logical, intent(in), optional :: south_halo
     logical, intent(in), optional :: north_halo
 
-    ! integer, allocatable :: buffer(:)
-    real, allocatable :: send_buffer1(:)
-    real, allocatable :: send_buffer2(:)
-    integer send_size
     logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
     integer t1, t2, t3, i, j, js, je, nx, mx, hx, hy, ierr
     integer send_req, recv_req
-    integer send_req1, recv_req1, send_req2, recv_req2, send_req3, recv_req3, send_req4, recv_req4
-    integer :: reqs1(1) 
-    integer :: reqs2(2)
-    integer :: reqs3(3)  
-    integer :: reqs4(4)
-    integer :: reqs5(5) 
-    integer :: reqs6(6)
-    integer :: reqs7(7)
-    integer :: reqs8(8)
-    integer counter
-    integer numreqs
-    integer iftest
-    integer :: status1(MPI_STATUS_SIZE*1)
-    integer :: status2(MPI_STATUS_SIZE*2)
-    integer :: status3(MPI_STATUS_SIZE*3)
-    integer :: status4(MPI_STATUS_SIZE*4)
-    integer :: status5(MPI_STATUS_SIZE*5)
-    integer :: status6(MPI_STATUS_SIZE*6)
-    integer :: status7(MPI_STATUS_SIZE*7)
-    integer :: status8(MPI_STATUS_SIZE*8)
     real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))
 
-
-
-    iftest = 1
     call perf_start('fill_halo_3d')
 
     west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
@@ -224,8 +195,6 @@ contains
     t3 = merge(1, 2, field%full_lev)
     hx = field%halo(1)%lon_hw
     hy = field%halo(1)%lat_hw
-    counter = 0
-    numreqs = 0
     if (field%full_lon) then
       nx = field%mesh%full_nlon
       mx = field%mesh%full_nlon / 2
@@ -241,17 +210,6 @@ contains
       je = field%mesh%half_jme
     end if
 
-    ! spherical_area_with_last_small_arc
-    ! swm 
-
-
-  ! omp {
-  !   #omp
-  !   for ()
-  ! }
-    ! swallow water mode
-
-  if (iftest .eq. 1) then 
     if (west_halo_opt) then
       call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
                         field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
@@ -292,313 +250,6 @@ contains
       call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
     end if
 
-
-  else 
-
-      if (west_halo_opt) then
-        call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
-                          field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
-                          proc%comm, MPI_STATUS_IGNORE, ierr)
-      end if
-
-      if (east_halo_opt) then
-        call MPI_SENDRECV(field%d, 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 32, &
-                          field%d, 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 32, &
-                          proc%comm, MPI_STATUS_IGNORE, ierr)  
-      end if
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    ! ! end if
-    !   ! deallocate(send_buffer1)
-    !   ! call MPI_Buffer_detach(recv_buffer, send_size + 1000 ,ierr)
-    ! end if
-
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    ! end if
-
-
-
-    call MPI_TYPE_SIZE(field%halo(north)%send_type_3d(t1,t2,t3), send_size, ierr)
-    allocate(send_buffer1(send_size * 2 + 1000))
-    ! ! allocate(send_buffer2(send_size + 1000))
-    call MPI_BUFFER_ATTACH(send_buffer1, send_size * 2 + 1000 ,ierr)
-
-      send_req3 = MPI_REQUEST_NULL; recv_req3 = MPI_REQUEST_NULL
-    if (south_halo_opt) then
-      ! send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_north_pole) then
-
-        call t_startf ('IBSEND')
-        call MPI_IbSEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-                       proc%comm, send_req3, ierr)
-        counter = counter + 1
-        call t_stopf ('IBSEND')
-      end if
-      if (.not. proc%at_south_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-                       proc%comm, recv_req3, ierr)
-        counter = counter + 1
-      end if
-      ! call MPI_WAIT(send_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-      ! deallocate(send_buffer1)
-      ! call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-    end if
-
-
-
-      send_req4 = MPI_REQUEST_NULL; recv_req4 = MPI_REQUEST_NULL
-    if (north_halo_opt) then
-      ! send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_south_pole) then
-        ! call MPI_BUFFER_ATTACH(send_buffer2, send_size + 1000 ,ierr)
-        call t_startf('IBSEND')
-        call MPI_IbSEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-                       proc%comm, send_req4, ierr)
-        counter = counter + 1
-        call t_stopf('IBSEND')
-      end if
-      if (.not. proc%at_north_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-                       proc%comm, recv_req4, ierr)
-        counter = counter + 1
-      end if
-
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req4, MPI_STATUS_IGNORE, ierr)
-    end if
-
-      ! call MPI_WAIT(send_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      if (counter .eq. 0) then 
-        counter = 0
-      else if (counter .eq. 1) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(1, reqs1, status1, ierr)
-      else if (counter .eq. 2) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(2, reqs2, status2, ierr)
-      else if (counter .eq. 3) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(3, reqs3, status3, ierr)
-      else if (counter .eq. 4) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(4, reqs4, status4, ierr)
-      end if
-
-      ! if (north_halo_opt .and. south_halo_opt) then 
-      !   if (.not. proc%at_south_pole) then 
-      !     counter = counter + 1
-      !     reqs4(counter) = send_req4
-      !     counter = counter + 1
-      !     reqs4(counter) = recv_req3
-      ! ! ! !   end if 
-      ! ! ! !   if (.not.  proc%at_north_pole) then 
-      !     counter = counter + 1
-      !     reqs4(counter) = send_req3
-      !     counter = counter + 1
-      !     reqs4(counter) = recv_req4
-      ! !   end if
-      !   MPI
-      ! else if (north_halo_opt) then 
-      ! if ((recv_req3 .ne. MPI_REQUEST_NULL) .and. (send_req3 .ne. MPI_REQUEST_NULL) .and. (recv_req4 .ne. MPI_REQUEST_NULL) .and. (send_req4 .ne. MPI_REQUEST_NULL)) then 
-      !   req3 += 
-      ! else if    
-      ! if (north)
-      ! else if (south_halo_opt) then 
-      ! if (north_halo_opt .and. north_halo_opt) then
-      ! if (.not. proc%at_north_pole ) then 
-
-      ! end if 
-      ! if (send_req3 == NULL ) 
-      ! call MPI_Waitany(4, reqs4, status4, ierr)
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req4, MPI_STATUS_IGNORE, ierr)
-
-      ! ! call MPI_WAIT(re)
-      call MPI_Buffer_detach(send_buffer1, send_size * 2 + 1000 ,ierr)
-      ! call MPI_Buffer_detach(send_buffer2, send_size + 1000 ,ierr)
-    deallocate(send_buffer1)
-    ! deallocate(send_buffer2)
-
-    ! call MPI_TYPE_SIZE(field%halo(north)%send_type_3d(t1,t2,t3), send_size, ierr)
-    ! allocate(send_buffer1(send_size + 1000))
-
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_BUFFER_ATTACH(send_buffer1, send_size + 1000 ,ierr)
-    !     call MPI_IbSEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-    !   end if
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-
-    !   ! deallocate(send_buffer1)
-    !   ! call MPI_Buffer_detach(recv_buffer, send_size + 1000 ,ierr)
-    ! end if
-
-
-    ! call MPI_TYPE_SIZE(field%halo(south)%send_type_3d(t1,t2,t3), send_size, ierr)
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   ! allocate(send_buffer2(send_size + 1000))
-
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_BUFFER_ATTACH(send_buffer1, send_size + 1000 ,ierr)
-    !     call MPI_IbSEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_Buffer_detach(send_buffer1, send_size + 1000, ierr)
-    !   end if
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-
-    ! end if
-
-    !   deallocate(send_buffer1)
-      
-
-
-
-
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (proc%at_south_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                   proc%comm, send_req, ierr)  
-    !     call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                   proc%comm, recv_req, ierr)
-    !     call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (.not. proc%at_north_pole .and. .not. proc%at_north_pole) then
-    !     call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                       field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                       proc%comm, MPI_STATUS_IGNORE, ierr)  
-    !   endif
-    ! end if
-
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (proc%at_north_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !     call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (.not. proc%at_north_pole .and. .not. proc%at_north_pole) then
-    !     call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                       field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                       proc%comm, MPI_STATUS_IGNORE, ierr)  
-    !   endif
-    ! end if
-    end if 
-
-
-
-
-    ! if (north_halo_opt .and. south_halo_opt) then
-
-
     if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
       call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
                         field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
@@ -640,699 +291,6 @@ contains
 
   end subroutine fill_halo_3d
 
-  subroutine fill_halo_start_3d(field, west_halo, east_halo, south_halo, north_halo, isstart, send_south_req, recv_south_req, send_north_req, recv_north_req)
-
-    type(latlon_field3d_type), intent(in) :: field
-    logical, intent(in), optional :: west_halo
-    logical, intent(in), optional :: east_halo
-    logical, intent(in), optional :: south_halo
-    logical, intent(in), optional :: north_halo
-    logical, intent(in) :: isstart
-    integer, intent(inout) :: send_south_req, recv_north_req, send_north_req, recv_south_req
-
-
-
-    ! integer, allocatable :: buffer(:)
-    real, allocatable :: send_buffer1(:)
-    real, allocatable :: send_buffer2(:)
-    integer send_size
-    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
-    integer t1, t2, t3, i, j, js, je, nx, mx, hx, hy, ierr
-    integer send_req, recv_req
-    integer send_req1, recv_req1, send_req2, recv_req2, send_req3, recv_req3, send_req4, recv_req4
-    integer :: reqs1(1) 
-    integer :: reqs2(2)
-    integer :: reqs3(3)  
-    integer :: reqs4(4)
-    integer :: reqs5(5) 
-    integer :: reqs6(6)
-    integer :: reqs7(7)
-    integer :: reqs8(8)
-    integer counter
-    integer numreqs
-    integer iftest
-    integer :: status1(MPI_STATUS_SIZE*1)
-    integer :: status2(MPI_STATUS_SIZE*2)
-    integer :: status3(MPI_STATUS_SIZE*3)
-    integer :: status4(MPI_STATUS_SIZE*4)
-    integer :: status5(MPI_STATUS_SIZE*5)
-    integer :: status6(MPI_STATUS_SIZE*6)
-    integer :: status7(MPI_STATUS_SIZE*7)
-    integer :: status8(MPI_STATUS_SIZE*8)
-    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))
-
-
-
-    iftest = 1
-    call perf_start('fill_halo_start_3d')
-    ! PRINT *, "Hello, Fortran!"
-
-    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
-    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
-    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
-    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo
-
-    t1 = merge(1, 2, field%full_lon)
-    t2 = merge(1, 2, field%full_lat)
-    t3 = merge(1, 2, field%full_lev)
-    hx = field%halo(1)%lon_hw
-    hy = field%halo(1)%lat_hw
-    counter = 0
-    numreqs = 0
-    if (field%full_lon) then
-      nx = field%mesh%full_nlon
-      mx = field%mesh%full_nlon / 2
-    else
-      nx = field%mesh%half_nlon
-      mx = field%mesh%half_nlon / 2
-    end if
-    if (field%full_lat) then
-      js = field%mesh%full_jms
-      je = field%mesh%full_jme
-    else
-      js = field%mesh%half_jms
-      je = field%mesh%half_jme
-    end if
-
-    ! spherical_area_with_last_small_arc
-    ! swm 
-
-
-  ! omp {
-  !   #omp
-  !   for ()
-  ! }
-    ! swallow water mode
-  send_south_req = MPI_REQUEST_NULL; recv_south_req  = MPI_REQUEST_NULL
-  send_north_req = MPI_REQUEST_NULL; recv_north_req  = MPI_REQUEST_NULL
-  if (iftest .eq. 1) then 
-    if (west_halo_opt) then
-      call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
-                        field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
-                        proc%comm, MPI_STATUS_IGNORE, ierr)
-    end if
-
-    if (east_halo_opt) then
-      call MPI_SENDRECV(field%d, 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 32, &
-                        field%d, 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 32, &
-                        proc%comm, MPI_STATUS_IGNORE, ierr)
-    end if
-
-    if (south_halo_opt) then
-      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_north_pole) then
-        call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-                       proc%comm, send_req, ierr)
-      end if
-      if (.not. proc%at_south_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-                       proc%comm, recv_req, ierr)
-      end if
-
-      send_south_req = send_req
-      recv_south_req = recv_req
-      ! call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    end if
-
-    if (north_halo_opt) then
-      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_south_pole) then
-        call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-                       proc%comm, send_req, ierr)
-      end if
-      if (.not. proc%at_north_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-                       proc%comm, recv_req, ierr)
-      end if
-      ! call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-      send_north_req = send_req
-      recv_north_req = recv_req
-      
-    end if
-
-
-  else 
-
-      if (west_halo_opt) then
-        call MPI_SENDRECV(field%d, 1, field%halo(east)%send_type_3d(t1,t2,t3), field%halo(east)%proc_id, 31, &
-                          field%d, 1, field%halo(west)%recv_type_3d(t1,t2,t3), field%halo(west)%proc_id, 31, &
-                          proc%comm, MPI_STATUS_IGNORE, ierr)
-      end if
-
-      if (east_halo_opt) then
-        call MPI_SENDRECV(field%d, 1, field%halo(west)%send_type_3d(t1,t2,t3), field%halo(west)%proc_id, 32, &
-                          field%d, 1, field%halo(east)%recv_type_3d(t1,t2,t3), field%halo(east)%proc_id, 32, &
-                          proc%comm, MPI_STATUS_IGNORE, ierr)  
-      end if
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    ! ! end if
-    !   ! deallocate(send_buffer1)
-    !   ! call MPI_Buffer_detach(recv_buffer, send_size + 1000 ,ierr)
-    ! end if
-
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    ! end if
-
-
-
-    call MPI_TYPE_SIZE(field%halo(north)%send_type_3d(t1,t2,t3), send_size, ierr)
-    allocate(send_buffer1(send_size * 2 + 1000))
-    ! ! allocate(send_buffer2(send_size + 1000))
-    call MPI_BUFFER_ATTACH(send_buffer1, send_size * 2 + 1000 ,ierr)
-
-      send_req3 = MPI_REQUEST_NULL; recv_req3 = MPI_REQUEST_NULL
-    if (south_halo_opt) then
-      ! send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_north_pole) then
-
-        call t_startf ('IBSEND')
-        call MPI_IbSEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-                       proc%comm, send_req3, ierr)
-        counter = counter + 1
-        call t_stopf ('IBSEND')
-      end if
-      if (.not. proc%at_south_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-                       proc%comm, recv_req3, ierr)
-        counter = counter + 1
-      end if
-      ! call MPI_WAIT(send_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-      ! deallocate(send_buffer1)
-      ! call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-    end if
-
-
-
-      send_req4 = MPI_REQUEST_NULL; recv_req4 = MPI_REQUEST_NULL
-    if (north_halo_opt) then
-      ! send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-      if (.not. proc%at_south_pole) then
-        ! call MPI_BUFFER_ATTACH(send_buffer2, send_size + 1000 ,ierr)
-        call t_startf('IBSEND')
-        call MPI_IbSEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-                       proc%comm, send_req4, ierr)
-        counter = counter + 1
-        call t_stopf('IBSEND')
-      end if
-      if (.not. proc%at_north_pole) then
-        call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-                       proc%comm, recv_req4, ierr)
-        counter = counter + 1
-      end if
-
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req4, MPI_STATUS_IGNORE, ierr)
-    end if
-
-      ! call MPI_WAIT(send_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      if (counter .eq. 0) then 
-        counter = 0
-      else if (counter .eq. 1) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs1(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(1, reqs1, status1, ierr)
-      else if (counter .eq. 2) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs2(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(2, reqs2, status2, ierr)
-      else if (counter .eq. 3) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs3(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(3, reqs3, status3, ierr)
-      else if (counter .eq. 4) then 
-          if (send_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = send_req4
-          end if
-          if (send_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = send_req3
-          end if
-          if (recv_req3 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = recv_req3
-          end if
-          if (recv_req4 .ne. MPI_REQUEST_NULL) then
-            numreqs = numreqs + 1
-            reqs4(numreqs) = recv_req4
-          end if
-        call MPI_Waitall(4, reqs4, status4, ierr)
-      end if
-
-      ! if (north_halo_opt .and. south_halo_opt) then 
-      !   if (.not. proc%at_south_pole) then 
-      !     counter = counter + 1
-      !     reqs4(counter) = send_req4
-      !     counter = counter + 1
-      !     reqs4(counter) = recv_req3
-      ! ! ! !   end if 
-      ! ! ! !   if (.not.  proc%at_north_pole) then 
-      !     counter = counter + 1
-      !     reqs4(counter) = send_req3
-      !     counter = counter + 1
-      !     reqs4(counter) = recv_req4
-      ! !   end if
-      !   MPI
-      ! else if (north_halo_opt) then 
-      ! if ((recv_req3 .ne. MPI_REQUEST_NULL) .and. (send_req3 .ne. MPI_REQUEST_NULL) .and. (recv_req4 .ne. MPI_REQUEST_NULL) .and. (send_req4 .ne. MPI_REQUEST_NULL)) then 
-      !   req3 += 
-      ! else if    
-      ! if (north)
-      ! else if (south_halo_opt) then 
-      ! if (north_halo_opt .and. north_halo_opt) then
-      ! if (.not. proc%at_north_pole ) then 
-
-      ! end if 
-      ! if (send_req3 == NULL ) 
-      ! call MPI_Waitany(4, reqs4, status4, ierr)
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req3, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(send_req4, MPI_STATUS_IGNORE, ierr)
-      ! call MPI_WAIT(recv_req4, MPI_STATUS_IGNORE, ierr)
-
-      ! ! call MPI_WAIT(re)
-      call MPI_Buffer_detach(send_buffer1, send_size * 2 + 1000 ,ierr)
-      ! call MPI_Buffer_detach(send_buffer2, send_size + 1000 ,ierr)
-    deallocate(send_buffer1)
-    ! deallocate(send_buffer2)
-
-    ! call MPI_TYPE_SIZE(field%halo(north)%send_type_3d(t1,t2,t3), send_size, ierr)
-    ! allocate(send_buffer1(send_size + 1000))
-
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_BUFFER_ATTACH(send_buffer1, send_size + 1000 ,ierr)
-    !     call MPI_IbSEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_Buffer_detach(send_buffer1, send_size + 1000 ,ierr)
-    !   end if
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                    proc%comm, send_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-
-    !   ! deallocate(send_buffer1)
-    !   ! call MPI_Buffer_detach(recv_buffer, send_size + 1000 ,ierr)
-    ! end if
-
-
-    ! call MPI_TYPE_SIZE(field%halo(south)%send_type_3d(t1,t2,t3), send_size, ierr)
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   ! allocate(send_buffer2(send_size + 1000))
-
-    !   if (.not. proc%at_south_pole) then
-    !     call MPI_BUFFER_ATTACH(send_buffer1, send_size + 1000 ,ierr)
-    !     call MPI_IbSEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_Buffer_detach(send_buffer1, send_size + 1000, ierr)
-    !   end if
-    !   if (.not. proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !   end if
-    !   call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-
-    ! end if
-
-    !   deallocate(send_buffer1)
-      
-
-
-
-
-
-    ! if (south_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (proc%at_south_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                   proc%comm, send_req, ierr)  
-    !     call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (proc%at_north_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                   proc%comm, recv_req, ierr)
-    !     call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (.not. proc%at_north_pole .and. .not. proc%at_north_pole) then
-    !     call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 33, &
-    !                       field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 33, &
-    !                       proc%comm, MPI_STATUS_IGNORE, ierr)  
-    !   endif
-    ! end if
-
-    ! if (north_halo_opt) then
-    !   send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
-    !   if (proc%at_north_pole) then
-    !     call MPI_ISEND(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                    proc%comm, send_req, ierr)
-    !     call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (proc%at_south_pole) then
-    !     call MPI_IRECV(field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                    proc%comm, recv_req, ierr)
-    !     call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
-    !   else if (.not. proc%at_north_pole .and. .not. proc%at_north_pole) then
-    !     call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 34, &
-    !                       field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 34, &
-    !                       proc%comm, MPI_STATUS_IGNORE, ierr)  
-    !   endif
-    ! end if
-    end if 
-
-
-
-
-    ! if (north_halo_opt .and. south_halo_opt) then
-
-
-    ! if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
-    !   call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
-    !                     field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
-    !                     proc%comm, MPI_STATUS_IGNORE, ierr)
-    !   ! Reverse array order.
-    !   tmp = field%d(:,js:0,:)
-    !   if (field%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
-    !     do j = js, 0
-    !       field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
-    !       field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,hy+js-j,:)
-    !     end do
-    !   else
-    !     do j = js, 0
-    !       field%d(:,j,:) = tmp(:,hy+js-j,:)
-    !     end do
-    !   end if
-    ! end if
-
-    ! if (north_halo_opt .and. proc%at_north_pole .and. field%halo_cross_pole) then
-    !   send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
-    !   call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
-    !                     field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
-    !                     proc%comm, MPI_STATUS_IGNORE, ierr)
-    !   ! Reverse array order.
-    !   tmp = field%d(:,je-hy+1:je,:)
-    !   if (field%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
-    !     do j = je - hy + 1, je
-    !       field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,je+1-j,:)
-    !       field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,je+1-j,:)
-    !     end do
-    !   else
-    !     do j = je - hy + 1, je
-    !       field%d(:,j,:) = tmp(:,je+1-j,:)
-    !     end do
-    !   end if
-    ! end if
-
-    call perf_stop('fill_halo_start_3d')
-
-  end subroutine fill_halo_start_3d
-
-
-  subroutine fill_halo_stop_3d(field, west_halo, east_halo, south_halo, north_halo, isstart, isstop, send_south_req, recv_south_req, send_north_req, recv_north_req)
-
-    type(latlon_field3d_type), intent(in) :: field
-    logical, intent(in), optional :: west_halo
-    logical, intent(in), optional :: east_halo
-    logical, intent(in), optional :: south_halo
-    logical, intent(in), optional :: north_halo
-    logical, intent(in) :: isstart, isstop
-    integer, intent(in) :: send_south_req, recv_north_req, send_north_req, recv_south_req
-
-
-    ! integer, allocatable :: buffer(:)
-    ! real, allocatable :: send_buffer1(:)
-    ! real, allocatable :: send_buffer2(:)
-    integer send_size
-    logical west_halo_opt, east_halo_opt, south_halo_opt, north_halo_opt
-    integer t1, t2, t3, i, j, js, je, nx, mx, hx, hy, ierr
-    integer send_req, recv_req
-    integer send_req1, recv_req1, send_req2, recv_req2, send_req3, recv_req3, send_req4, recv_req4
-    integer :: reqs1(1) 
-    integer :: reqs2(2)
-    integer :: reqs3(3)  
-    integer :: reqs4(4)
-    integer :: reqs5(5) 
-    integer :: reqs6(6)
-    integer :: reqs7(7)
-    integer :: reqs8(8)
-    integer counter
-    integer numreqs
-    integer iftest
-    integer :: status1(MPI_STATUS_SIZE*1)
-    integer :: status2(MPI_STATUS_SIZE*2)
-    integer :: status3(MPI_STATUS_SIZE*3)
-    integer :: status4(MPI_STATUS_SIZE*4)
-    integer :: status5(MPI_STATUS_SIZE*5)
-    integer :: status6(MPI_STATUS_SIZE*6)
-    integer :: status7(MPI_STATUS_SIZE*7)
-    integer :: status8(MPI_STATUS_SIZE*8)
-    real(r8) tmp(size(field%d,1),field%halo(1)%lat_hw,size(field%d,3))
-
-
-
-    iftest = 1
-    call perf_start('fill_halo_stop_3d')
-
-    west_halo_opt  = .true. ; if (present(west_halo )) west_halo_opt  = west_halo
-    east_halo_opt  = .true. ; if (present(east_halo )) east_halo_opt  = east_halo
-    south_halo_opt = .true. ; if (present(south_halo)) south_halo_opt = south_halo
-    north_halo_opt = .true. ; if (present(north_halo)) north_halo_opt = north_halo
-
-    t1 = merge(1, 2, field%full_lon)
-    t2 = merge(1, 2, field%full_lat)
-    t3 = merge(1, 2, field%full_lev)
-    hx = field%halo(1)%lon_hw
-    hy = field%halo(1)%lat_hw
-    counter = 0
-    numreqs = 0
-    if (field%full_lon) then
-      nx = field%mesh%full_nlon
-      mx = field%mesh%full_nlon / 2
-    else
-      nx = field%mesh%half_nlon
-      mx = field%mesh%half_nlon / 2
-    end if
-    if (field%full_lat) then
-      js = field%mesh%full_jms
-      je = field%mesh%full_jme
-    else
-      js = field%mesh%half_jms
-      je = field%mesh%half_jme
-    end if
-
-    counter = 0
-
-    if (south_halo_opt) then
-      if (.not. proc%at_north_pole) then
-        counter = counter + 1
-      end if
-      if (.not. proc%at_south_pole) then
-        counter = counter + 1
-      end if
-    end if
-
-    if (north_halo_opt) then
-      if (.not. proc%at_south_pole) then
-        counter = counter + 1
-      end if
-      if (.not. proc%at_north_pole) then
-        counter = counter + 1
-      end if
-    end if
-
-
-    if (counter .eq. 0) then 
-      counter = 0
-    else if (counter .eq. 1) then 
-        if (send_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs1(numreqs) = send_south_req
-        end if
-        if (recv_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs1(numreqs) = recv_south_req
-        end if
-        if (send_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs1(numreqs) = send_north_req
-        end if
-        if (recv_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs1(numreqs) = recv_north_req
-        end if
-      call MPI_Waitall(1, reqs1, status1, ierr)
-    else if (counter .eq. 2) then 
-        if (send_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs2(numreqs) = send_south_req
-        end if
-        if (recv_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs2(numreqs) = recv_south_req
-        end if
-        if (send_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs2(numreqs) = send_north_req
-        end if
-        if (recv_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs2(numreqs) = recv_north_req
-        end if
-      call MPI_Waitall(2, reqs2, status2, ierr)
-    else if (counter .eq. 3) then 
-        if (send_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs3(numreqs) = send_south_req
-        end if
-        if (recv_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs3(numreqs) = recv_south_req
-        end if
-        if (send_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs3(numreqs) = send_north_req
-        end if
-        if (recv_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs3(numreqs) = recv_north_req
-        end if
-      call MPI_Waitall(3, reqs3, status3, ierr)
-    else if (counter .eq. 4) then 
-        if (send_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs4(numreqs) = send_south_req
-        end if
-        if (recv_south_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs4(numreqs) = recv_south_req
-        end if
-        if (send_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs4(numreqs) = send_north_req
-        end if
-        if (recv_north_req .ne. MPI_REQUEST_NULL) then
-          numreqs = numreqs + 1
-          reqs4(numreqs) = recv_north_req
-        end if
-      call MPI_Waitall(4, reqs4, status4, ierr)
-    end if
-
-
-    if (south_halo_opt .and. proc%at_south_pole .and. field%halo_cross_pole) then
-      call MPI_SENDRECV(field%d, 1, field%halo(south)%send_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
-                        field%d, 1, field%halo(south)%recv_type_3d(t1,t2,t3), field%halo(south)%proc_id, 35, &
-                        proc%comm, MPI_STATUS_IGNORE, ierr)
-      ! Reverse array order.
-      tmp = field%d(:,js:0,:)
-      if (field%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
-        do j = js, 0
-          field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,hy+js-j,:)
-          field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,hy+js-j,:)
-        end do
-      else
-        do j = js, 0
-          field%d(:,j,:) = tmp(:,hy+js-j,:)
-        end do
-      end if
-    end if
-
-    if (north_halo_opt .and. proc%at_north_pole .and. field%halo_cross_pole) then
-      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
-      call MPI_SENDRECV(field%d, 1, field%halo(north)%send_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
-                        field%d, 1, field%halo(north)%recv_type_3d(t1,t2,t3), field%halo(north)%proc_id, 36, &
-                        proc%comm, MPI_STATUS_IGNORE, ierr)
-      ! Reverse array order.
-      tmp = field%d(:,je-hy+1:je,:)
-      if (field%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
-        do j = je - hy + 1, je
-          field%d(   1:mx,j,:) = tmp(hx+1+mx:hx+nx,je+1-j,:)
-          field%d(mx+1:nx,j,:) = tmp(hx+1   :hx+mx,je+1-j,:)
-        end do
-      else
-        do j = je - hy + 1, je
-          field%d(:,j,:) = tmp(:,je+1-j,:)
-        end do
-      end if
-    end if
-
-    call perf_stop('fill_halo_stop_3d')
-
-  end subroutine fill_halo_stop_3d
-
-
   subroutine fill_halo_4d(field, i4, west_halo, east_halo, south_halo, north_halo)
 
     type(latlon_field4d_type), intent(in) :: field
```
**gmcore/src/physics/physics_mod.F90**
```fortran
diff --git a/gmcore/src/physics/physics_mod.F90 b/gmcore/src/physics/physics_mod.F90
index dc6da4b..4311e4b 100644
--- a/gmcore/src/physics/physics_mod.F90
+++ b/gmcore/src/physics/physics_mod.F90
@@ -19,7 +19,6 @@ module physics_mod
   use dp_coupling_mod
   use latlon_parallel_mod
   use simple_physics_driver_mod
-  use perf_mod
 #ifdef HAS_CAM
   use cam_physics_driver_mod
 #endif
@@ -220,7 +219,7 @@ contains
     real(r8), intent(in) :: dt
 
     integer i, j, k
-    call t_startf ( 'physics_update_dynamics' )
+
     associate (mesh  => block%mesh               , &
                dudt  => block%aux%dudt_phys      , & ! in
                dvdt  => block%aux%dvdt_phys      , & ! in
@@ -255,7 +254,6 @@ contains
     end do
     call fill_halo(pt)
     end associate
-    call t_stopf ( 'physics_update_dynamics' )
 
   end subroutine physics_update_dynamics
```
**gmcore/src/utils/perf_mod.F90**
```fortran
diff --git a/gmcore/src/utils/perf_mod.F90 b/gmcore/src/utils/perf_mod.F90
index a9ef18c..0deb00d 100644
--- a/gmcore/src/utils/perf_mod.F90
+++ b/gmcore/src/utils/perf_mod.F90
@@ -44,15 +44,7 @@ contains
     end if
     ierr = gptlstart('total')
 #endif
-    ! write(*, '(A)') 'hello wwwww'
-    ! ierr = gptlsetoption(gptloverhead, 0)
-    ! ierr = gptlsetoption(gptlpercent, 0)
-    ! ierr = gptlsetoption(gptlabort_on_error, 1)
-    ! ierr = gptlsetoption(gptlsync_mpi, 1)
-    ! if (gptlinitialize() < 0) then
-    !   stop 'Failed to initialize GPTL!'
-    ! end if
-    ! ierr = gptlstart('total')
+
     initialized = .true.
 
   end subroutine perf_init
@@ -64,7 +56,6 @@ contains
     integer ierr
 
 #ifdef HAS_GPTL
-    ! write(*, '(A)') 'has gptl'
     ierr = gptlstart(name)
 #endif
 
@@ -93,12 +84,6 @@ contains
       stop 'Failed to call finalize GPTL!'
     end if
 #endif
-    ! write(*, '(A)') 'hello genshin'
-    ! ierr = gptlstop('total')
-    ! ierr = gptlpr(0)
-    ! if (gptlfinalize() /= 0) then
-    !   stop 'Failed to call finalize GPTL!'
-    ! end if
 
   end subroutine perf_final
```
**gmcore/src/utils/time_mod.F90**
```fortran
diff --git a/gmcore/src/utils/time_mod.F90 b/gmcore/src/utils/time_mod.F90
index a42af35..eb777b5 100644
--- a/gmcore/src/utils/time_mod.F90
+++ b/gmcore/src/utils/time_mod.F90
@@ -5,7 +5,6 @@ module time_mod
   use container
   use flogger
   use const_mod
-  use perf_mod
   use namelist_mod, start_time_array => start_time, end_time_array => end_time
 
   implicit none
@@ -157,7 +156,6 @@ contains
 
     type(hash_table_iterator_type) iter
 
-    call t_startf ( 'time_advance' )
     ! Update alerts.
     iter = hash_table_iterator(alerts)
     do while (.not. iter%ended())
@@ -183,8 +181,6 @@ contains
     end if
     curr_time_str = curr_time%isoformat()
 
-    call t_stopf ( 'time_advance' )
-
   end subroutine time_advance
 
   subroutine time_fast_forward(time_value, time_units)
```