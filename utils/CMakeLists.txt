set(STELLA_UTILS_SOURCES_f90
  convert.f90
  fft_work.f90
  gauss_quad.f90
  linear_solve.f90
  mt19937.f90
  redistribute.f90
  smooth_step.f90
  sort.f90
  spl.f90
  text_options.f90
  )

set(STELLA_UTILS_SOURCES_fpp
  command_line.fpp
  constants.fpp
  file_utils.fpp
  job_manage.fpp
  mp.fpp
  mp_lu_decomposition.fpp
  ran.fpp
  spfunc.fpp
  system_fortran.fpp
  )
set_source_files_properties(
  ${STELLA_UTILS_SOURCES_fpp}
  PROPERTIES
  Fortran_PREPROCESS ON
  )

add_library(stella_utils
  ${STELLA_UTILS_SOURCES_f90}
  ${STELLA_UTILS_SOURCES_fpp}
  )

# Location of built and installed .mod files
set(stella_utils_module_directory ${CMAKE_CURRENT_BINARY_DIR}/mod)
set(stella_utils_module_install_directory "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING
  "Fortran module installation directory")

# Make sure all source files are compiled in free-format
set_target_properties(stella_utils PROPERTIES
  Fortran_FORMAT FREE
  Fortran_MODULE_DIRECTORY ${stella_utils_module_directory}
  )

target_include_directories(stella_utils PUBLIC
  $<BUILD_INTERFACE:${stella_utils_module_directory}>
  $<INSTALL_INTERFACE:${stella_utils_module_install_directory}>
  )

if (STELLA_ENABLE_F200X AND STELLA_ENABLE_LOCAL_SPFUNC)
  message(FATAL_ERROR "STELLA_ENABLE_F200X is incompatible with STELLA_ENABLE_LOCAL_SPFUNC")
endif()

target_link_libraries(stella_utils PUBLIC ${_fftw_lib})
target_compile_definitions(stella_utils PUBLIC FFT=_FFTW3_)

if (STELLA_ENABLE_NETCDF)
  target_sources(stella_utils PRIVATE netcdf_utils.fpp)
  set_source_files_properties(netcdf_utils.fpp PROPERTIES
    Fortran_PREPROCESS ON
    Fortran_FORMAT FREE
    )
  find_package(netCDFFortran REQUIRED)
  target_link_libraries(stella_utils PUBLIC netCDF::netcdff)
  target_compile_definitions(stella_utils PRIVATE NETCDF=_DEFAULT_)
endif()

if (STELLA_ENABLE_DOUBLE)
  target_compile_options(stella_utils PUBLIC
    $<$<Fortran_COMPILER_ID:GNU>:
      -fdefault-real-8 -fdefault-double-8 >
    $<$<Fortran_COMPILER_ID:Intel>:
      -r8 >
    )
  set(STELLA_DEFAULT_REAL_KIND "double")
else()
  set(STELLA_DEFAULT_REAL_KIND "single")
  target_compile_definitions(stella_utils PRIVATE SINGLE_PRECISION)
endif()

# Type of special mathematical functions
if (STELLA_ENABLE_LOCAL_SPFUNC)
  set(SPFUNC _SPLOCAL_)
  set(STELLA_SPFUNC_TYPE "Local (stella_utils backport implementations)")
elseif(STELLA_ENABLE_NAGLIB)
  set(SPFUNC _SPNAG_)
  set(STELLA_SPFUNC_TYPE "NAG")
else()
  set(SPFUNC _SPF200X_)
  set(STELLA_SPFUNC_TYPE "F2003/F2008 intrinsics")
endif()
target_compile_definitions(stella_utils PRIVATE SPFUNC=${SPFUNC})

if (STELLA_ENABLE_NAGLIB)
  if (STELLA_ENABLE_DOUBLE)
    target_compile_definitions(stella_utils PRIVATE NAG_PREC=_NAGDBLE_)
  else()
    target_compile_definitions(stella_utils PRIVATE NAG_PREC=_NAGSNGL_)
  endif()
endif()

if (STELLA_ENABLE_POSIX)
  target_compile_definitions(stella_utils PRIVATE POSIX)
endif()

if (STELLA_ENABLE_F200X)
  target_compile_definitions(stella_utils PRIVATE F200X_INTRINSICS)
endif()

if (STELLA_ENABLE_RANDOM_MT)
  target_compile_definitions(stella_utils PRIVATE RANDOM=_RANMT_)
endif()

# Ensure ANSI_CPP defined if using intel (required for preprocessing)
target_compile_definitions(stella_utils PRIVATE
  $<$<Fortran_COMPILER_ID:Intel>:ANSI_CPP>)

# Ensure FCOMPILER is set
target_compile_definitions(stella_utils PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:FCOMPILER=_GFORTRAN_>
  $<$<Fortran_COMPILER_ID:Intel>:FCOMPILER=_INTEL_>
  )

find_package(MPI)
target_link_libraries(stella_utils PUBLIC MPI::MPI_Fortran)
target_compile_definitions(stella_utils PRIVATE MPI)
