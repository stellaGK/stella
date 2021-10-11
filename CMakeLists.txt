cmake_minimum_required(VERSION 3.16)
project(stella
  VERSION 0.3
  LANGUAGES Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Sets default build type to RelWithDebInfo
include(BuildType)

############################################################
# Build options
option(STELLA_ENABLE_LAPACK "Enable LAPACK wrappers" ON)
option(STELLA_ENABLE_FFT "Enable FFTs" ON)
option(STELLA_ENABLE_NETCDF "Enable NetCDF" ON)
option(STELLA_ENABLE_DOUBLE "Promotes precisions of real and complex to double" ON)
option(STELLA_ENABLE_PETSC "Enable PETSc" OFF)

# Options for utils internal library
option(STELLA_ENABLE_LOCAL_SPFUNC "Enable local special functions" OFF)
option(STELLA_ENABLE_NAGLIB "Use the NAG library" OFF)
option(STELLA_ENABLE_POSIX "Enable POSIX functions for command line functionality" OFF)
option(STELLA_ENABLE_RANDOM_MT "Enable Mersenne-Twister rng (not recommended)" OFF)
option(STELLA_ENABLE_F200X "Enable use of F2003/F2008 functionality" ON)

# Options for testing
option(STELLA_ENABLE_TESTS "Enable pFUnit-based tests" ON)
option(STELLA_DOWNLOAD_PFUNIT "Download and compile pFUnit at configure time" ON)

############################################################
# Sources

# Sources that don't need preprocessing
set(STELLA_SOURCES_f90
  common_types.f90
  define.inc
  dissipation.f90
  dist_fn.f90
  dist_fn_arrays.f90
  dist_redistribute.f90
  euterpe_interface.f90
  extended_zgrid.f90
  fields_arrays.f90
  finite_differences.f90
  flow_shear.f90
  g_tofrom_h.f90
  gyro_averages.f90
  init_g.f90
  kt_grids.f90
  mirror_terms.f90
  neoclassical_terms.f90
  parallel_streaming.f90
  physics_flags.f90
  physics_parameters.f90
  run_parameters.f90
  species.f90
  stella_diagnostics.f90
  stella_layouts.f90
  stella_time.f90
  stella_transforms.f90
  time_advance.f90
  volume_averages.f90
  vpamu_grids.f90
  zgrid.f90
  geo/inputprofiles_interface.f90
  geo/millerlocal.f90
  geo/stella_geometry.f90
  geo/vmec_geo.f90
  geo/vmec_interface/vmec_to_stella_geometry_interface.f90
  geo/vmec_interface/fzero.f90
)

# Sources that _do_ need preprocessing
set(STELLA_SOURCES_fpp
  fields.fpp
  multibox.fpp
  response_matrix.fpp
  sfincs_interface.fpp
  sources.fpp
  stella_io.fpp
  stella_save.fpp
)
set_source_files_properties(${STELLA_SOURCES_fpp} PROPERTIES
  Fortran_PREPROCESS ON
  )

# Internal library that we can link tests against
add_library(libstella OBJECT
  ${STELLA_SOURCES_f90}
  ${STELLA_SOURCES_fpp}
  )

# Location of built and installed .mod files
set(libstella_module_directory ${CMAKE_CURRENT_BINARY_DIR}/mod)
set(libstella_module_install_directory "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING
  "Fortran module installation directory")

# Make sure all source files are compiled in free-format
set_target_properties(libstella PROPERTIES
  Fortran_FORMAT FREE
  Fortran_MODULE_DIRECTORY ${libstella_module_directory}
  )

target_include_directories(libstella PUBLIC
  $<BUILD_INTERFACE:${libstella_module_directory}>
  $<INSTALL_INTERFACE:${libstella_module_install_directory}>
  )

# Don't worry about line length
target_compile_options(libstella PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>
  )

# Main executable
add_executable(stella
  stella.f90
  )
target_link_libraries(stella PRIVATE libstella)

############################################################
# Dependencies and build options

# This might not be entirely sensible, but helps CMake to find the
# correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
find_package(MPI REQUIRED)
target_link_libraries(libstella PUBLIC MPI::MPI_Fortran)
target_compile_definitions(libstella PUBLIC MPI)

if (STELLA_ENABLE_LAPACK)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(LAPACK REQUIRED)
    target_link_libraries(libstella PUBLIC "${LAPACK_LIBRARIES}")
  endif()
  target_compile_definitions(libstella PUBLIC LAPACK)
endif()

if (STELLA_ENABLE_FFT)
  if (STELLA_ENABLE_DOUBLE)
    set(_fftw_component DOUBLE_LIB)
    set(_fftw_lib FFTW::Double)
  else()
    set(_fftw_component FLOAT_LIB)
    set(_fftw_lib FFTW::Float)
  endif()
  find_package(FFTW REQUIRED COMPONENTS ${_fftw_component})
  target_link_libraries(libstella PUBLIC ${_fftw_lib})
  target_compile_definitions(libstella PUBLIC FFT=_FFTW3_)
else()
  message(WARNING "STELLA_ENABLE_FFT is off -- Be sure that nonlinear run makes no sense")
endif()

if (STELLA_ENABLE_NETCDF)
  find_package(netCDFFortran REQUIRED)
  target_link_libraries(libstella PUBLIC netCDF::netcdff)
  target_compile_definitions(libstella PUBLIC NETCDF=_DEFAULT_)
endif()

if (STELLA_ENABLE_DOUBLE)
  target_compile_options(libstella PUBLIC
    $<$<Fortran_COMPILER_ID:GNU>:
      -fdefault-real-8 -fdefault-double-8 >
    $<$<Fortran_COMPILER_ID:Intel>:
      -r8 >
    $<$<Fortran_COMPILER_ID:Cray>:
      -s real64 >
    )
  set(STELLA_DEFAULT_REAL_KIND "double")
else()
  set(STELLA_DEFAULT_REAL_KIND "single")
  target_compile_definitions(libstella PUBLIC SINGLE_PRECISION)
endif()

if (STELLA_ENABLE_F200X)
  target_compile_definitions(libstella PUBLIC ISO_C_BINDING)
endif()

if (STELLA_ENABLE_PETSC)
  find_package(PETSc REQUIRED)
  target_link_libraries(libstella PUBLIC PETSc::PETSc)
endif()

# Ensure FCOMPILER is set
target_compile_definitions(libstella PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:FCOMPILER=_GFORTRAN_>
  $<$<Fortran_COMPILER_ID:Intel>:FCOMPILER=_INTEL_>
  $<$<Fortran_COMPILER_ID:Cray>:FCOMPILER=_CRAY_>
  )

############################################################
# Add some extra flags for Debug configs

target_compile_options(libstella PRIVATE
  $<$<CONFIG:Debug>:
    $<$<Fortran_COMPILER_ID:GNU>:
      -g -Wall -fimplicit-none -fbounds-check >
    $<$<Fortran_COMPILER_ID:Intel>:
      -g -implicitnone -warn all -nogen-interfaces -CB -traceback >
    $<$<Fortran_COMPILER_ID:Cray>:
      -g -Rb >
    >
  )

##################################################
# Internal libraries

add_subdirectory(geo/vmec_interface)
add_subdirectory(utils)
target_link_libraries(libstella PUBLIC mini_libstell stella_utils)

############################################################
# Tests

if (STELLA_ENABLE_TESTS)
  if (STELLA_DOWNLOAD_PFUNIT)
    message(STATUS "Downloading pFUnit...")
    include(FetchContent)
    set(STELLA_DOWNLOAD_PFUNIT_URL "https://github.com/Goddard-Fortran-Ecosystem/pFUnit"
      CACHE STRING "URL of pFUnit repository")

    set(STELLA_PFUNIT_LOCATION "${STELLA_DOWNLOAD_PFUNIT_URL}")

    FetchContent_Declare(
      pfunit
      GIT_REPOSITORY ${STELLA_DOWNLOAD_PFUNIT_URL}
      GIT_TAG        main
      )

    set(SKIP_OPENMP YES CACHE BOOL "" FORCE)
    # Workaround for bug in current main at time of commit
    set(BUILD_SHARED ${BUILD_SHARED_LIBS})
    if (NOT CMAKE_Fortran_FLAGS)
      if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
        # Workaround for an issue in pFUnit: preprocessor macro
        # expansion of the current source directory, which may
        # truncate the line
        set(CMAKE_Fortran_FLAGS "-ffree-line-length-none" CACHE STRING "" FORCE)
        message(STATUS "Setting CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
      endif()
    endif()
    FetchContent_MakeAvailable(pfunit)
  else()
    find_package(PFUNIT REQUIRED)
    set(STELLA_PFUNIT_LOCATION "${PFUNIT_ROOT}")
  endif()
  include(CTest)
  enable_testing()

  include(GKfunctions)

  set(STELLA_CTEST_CASES)
  stella_add_standard_tests(TEST_NAMES
    simple
    utils_convert
    )
  add_custom_target(check
    COMMAND ctest --output-on-failure
    DEPENDS ${STELLA_CTEST_CASES}
    )
else()
  set(STELLA_PFUNIT_LOCATION "Not used")
endif()

##################################################
# Configure summary

# The optimisation and/or debug flags are in the CMAKE_Fortran_FLAGS*
# variables. We need both the common flags as well as those for the
# build type actually being used. Note: this might behave weirdly on
# Windows. Might need to expand CMAKE_CONFIGURATION_TYPES instead?
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE_UPPER "DEBUG")
else()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
endif()
string(JOIN " " STELLA_COMPILE_BUILD_FLAGS
  "${CMAKE_Fortran_FLAGS}"
  "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}"
  )

message("
   ---------------------------------
     stella Configuration Summary
   ---------------------------------

   LAPACK support           : ${STELLA_ENABLE_LAPACK}
   FFTW support             : ${STELLA_ENABLE_FFT}
   NetCDF support           : ${STELLA_ENABLE_NETCDF}
   Default real kind        : ${STELLA_DEFAULT_REAL_KIND}
   POSIX support            : ${STELLA_ENABLE_POSIX}
   Use F2003/8 intrinsics   : ${STELLA_ENABLE_F200X}
   Compilation flags        : ${STELLA_COMPILE_BUILD_FLAGS}
   Build type               : ${CMAKE_BUILD_TYPE}
   Compiler                 : ${CMAKE_Fortran_COMPILER}
")
