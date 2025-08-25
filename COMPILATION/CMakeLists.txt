
############################################################
#  Create the internal library <libstella>
############################################################

cmake_minimum_required(VERSION 3.16)
project(stella
  VERSION 0.3
  LANGUAGES Fortran)

# Cmake files
set(CMAKE_MODULE_PATH
"${CMAKE_SOURCE_DIR}/COMPILATION/cmake" 
${CMAKE_MODULE_PATH})

# Sets default build type to RelWithDebInfo
include(BuildType)

############################################################
#  Create some options for the compilation
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
option(STELLA_ENABLE_TESTS "Enable pFUnit-based tests" OFF)
option(STELLA_DOWNLOAD_PFUNIT "Download and compile pFUnit at configure time" ON)
option(STELLA_UPDATE_GIT_SUBMODULES "Check submodules are up-to-date during build" ON)

############################################################
# Define the source files for stella
############################################################

# Sources that don't need preprocessing
set(STELLA_SOURCES_f90 
  ../STELLA_CODE/utils/define.inc
  ../STELLA_CODE/calculations/calculations_kxky.f90
  ../STELLA_CODE/calculations/calculations_checksum.f90
  ../STELLA_CODE/calculations/calculations_kxky_derivatives.f90
  ../STELLA_CODE/calculations/calculations_redistribute.f90
  ../STELLA_CODE/calculations/calculations_finite_differences.f90
  ../STELLA_CODE/calculations/calculations_tofrom_ghf.f90
  ../STELLA_CODE/calculations/calculations_gyro_averages.f90
  ../STELLA_CODE/calculations/calculations_transforms.f90
  ../STELLA_CODE/calculations/calculations_volume_averages.f90
  ../STELLA_CODE/diagnostics/diagnostics.f90
  ../STELLA_CODE/diagnostics/diagnostics_distribution.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_fluxtube.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_fullfluxsurface.f90 
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_radialvariation.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes.f90 
  ../STELLA_CODE/diagnostics/diagnostics_moments.f90
  ../STELLA_CODE/diagnostics/diagnostics_omega.f90
  ../STELLA_CODE/diagnostics/diagnostics_potential.f90
  ../STELLA_CODE/dissipation/dissipation_coll_dougherty.f90
  ../STELLA_CODE/dissipation/dissipation_coll_fokkerplanck.f90
  ../STELLA_CODE/dissipation/dissipation.f90
  ../STELLA_CODE/dissipation/dissipation_hyper.f90
  ../STELLA_CODE/fields/initialise_distribution_fn.f90
  ../STELLA_CODE/geometry/geometry.f90
  ../STELLA_CODE/geometry/geometry_miller.f90
  ../STELLA_CODE/geometry/geometry_vmec.f90
  ../STELLA_CODE/geometry/geometry_vmec_read_netCDF_file.f90 
  ../STELLA_CODE/geometry/geometry_inputprofiles_interface.f90
  ../STELLA_CODE/geometry/geometry_zpinch.f90
  ../STELLA_CODE/grids/arrays_constants.f90
  ../STELLA_CODE/grids/arrays_distribution_fn.f90
  ../STELLA_CODE/grids/arrays_gyro_averages.f90
  ../STELLA_CODE/grids/common_types.f90
  ../STELLA_CODE/grids/euterpe_interface.f90
  ../STELLA_CODE/grids/grids_extended_zgrid.f90
  ../STELLA_CODE/grids/grids_kxky.f90
  ../STELLA_CODE/grids/species.f90
  ../STELLA_CODE/grids/store_arrays_distribution_fn.f90
  ../STELLA_CODE/grids/store_arrays_fields.f90
  ../STELLA_CODE/grids/store_arrays_useful.f90 
  ../STELLA_CODE/grids/velocity_grids.f90
  ../STELLA_CODE/grids/write_radial_grid.f90
  ../STELLA_CODE/grids/z_grid.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_ffs_solve.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_radial_variation.f90
  ../STELLA_CODE/gyrokinetic_equation/reset_timestep.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_flow_shear.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_implicit_solve.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_mirror.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_parallel_streaming.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_nonlinearity.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_time_advance.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_drifts.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_fields.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_dissipation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_physics_parameters.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_numerical_parameters.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_kxky_grid.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_z_grid.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_species.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_velocity_grids.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_geometry.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_diagnostics.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_debug.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_parallelisation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_sources.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_neoclassical_input.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_radially_global.f90
  ../STELLA_CODE/neoclassical/euterpe_interface.f90
  ../STELLA_CODE/neoclassical/neoclassical_terms.f90
  ../STELLA_CODE/parameters/debug_flags.f90
  ../STELLA_CODE/parameters/parameters_diagnostics.f90
  ../STELLA_CODE/parameters/kxky_grid_parameters.f90
  ../STELLA_CODE/parameters/parameters_multibox.f90
  ../STELLA_CODE/parameters/numerical_parameters.f90
  ../STELLA_CODE/parameters/physics_parameters.f90
  ../STELLA_CODE/radial_variation/multibox.f90
  ../STELLA_CODE/stella_calculations/stella_layouts.f90
  ../STELLA_CODE/stella_calculations/stella_time.f90
)

# Sources that _do_ need preprocessing (because they have MPI flags, ...)
set(STELLA_SOURCES_fpp
  ../STELLA_CODE/stella_calculations/stella_io.fpp
  ../STELLA_CODE/stella_calculations/stella_save.fpp
  ../STELLA_CODE/fields/fields.fpp
  ../STELLA_CODE/fields/fields_fluxtube.fpp
  ../STELLA_CODE/fields/fields_electromagnetic.fpp
  ../STELLA_CODE/fields/fields_fullfluxsurface.fpp
  ../STELLA_CODE/fields/fields_collisions.fpp
  ../STELLA_CODE/fields/fields_radialvariation.fpp
  ../STELLA_CODE/gyrokinetic_equation/response_matrix.fpp
  ../STELLA_CODE/gyrokinetic_equation/sources.fpp
  ../STELLA_CODE/neoclassical/sfincs_interface.fpp
)

set_source_files_properties(${STELLA_SOURCES_fpp} PROPERTIES
  Fortran_PREPROCESS ON
  )

# Internal library <libstella> that we can link tests against
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

# Specifies include directories to use when compiling a given target
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
  ../STELLA_CODE/stella.f90
  )

# Use <libstella> when linking the <stella> target
target_link_libraries(stella PRIVATE libstella)

# The include.inc file is located in EXTERNALS/utils/
include_directories("${CMAKE_SOURCE_DIR}/STELLA_CODE/utils")

# When target `stella` is built, a symlink will be created
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_custom_command(TARGET stella POST_BUILD
    COMMAND echo Create symbolic link to stella in the directory above COMPILATION/build_cmake
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/stella ${CMAKE_CURRENT_BINARY_DIR}/../../stella
    DEPENDS stella
    VERBATIM ON
)
add_custom_command(TARGET stella POST_BUILD
    COMMAND echo Create symbolic link to stella in the directory above COMPILATION
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/stella ${CMAKE_CURRENT_BINARY_DIR}/../../../stella
    DEPENDS stella
    VERBATIM ON
)

############################################################
# Dependencies and build options
############################################################

# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
# Update submodules as needed
if(STELLA_UPDATE_GIT_SUBMODULES)
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

# This might not be entirely sensible, but helps CMake to find the
# correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
find_package(MPI REQUIRED)
target_link_libraries(libstella PUBLIC MPI::MPI_Fortran)
target_compile_definitions(libstella PUBLIC MPI)

############################################################
# Toggle the options defined at the start of this script
############################################################

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
  add_subdirectory(../EXTERNALS/neasyf ${CMAKE_CURRENT_BINARY_DIR}/EXTERNALS/neasyf)
  target_link_libraries(libstella PUBLIC netCDF::netcdff neasyf::neasyf)
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

# Link the <git_version> program
add_subdirectory(../EXTERNALS/git_version ${CMAKE_CURRENT_BINARY_DIR}/EXTERNALS/git_version)
target_link_libraries(libstella PUBLIC fortran_git::fortran_git)

# Ensure FCOMPILER is set
target_compile_definitions(libstella PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:FCOMPILER=_GFORTRAN_>
  $<$<Fortran_COMPILER_ID:Intel>:FCOMPILER=_INTEL_>
  $<$<Fortran_COMPILER_ID:Cray>:FCOMPILER=_CRAY_>
  )
  

# Find hdf5 packages (their CMake file needs the C language)
# These lines are only needed for Mac systems on Github since January 2025
# For some reason CMake could not find the hdf5 libraries anymore after a Github update
ENABLE_LANGUAGE(C)
find_package(HDF5 REQUIRED)
target_link_libraries(libstella PUBLIC ${HDF5_LIBRARIES}) 

############################################################
# Add some extra flags for Debug configs
############################################################

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

############################################################
# Internal libraries: <stella_utils> and <mini_libstell>
############################################################

add_subdirectory(../STELLA_CODE/utils ${CMAKE_CURRENT_BINARY_DIR}/STELLA_CODE/utils)
add_subdirectory(../EXTERNALS/mini_libstell ${CMAKE_CURRENT_BINARY_DIR}/EXTERNALS/mini_libstell)
target_link_libraries(libstella PUBLIC mini_libstell stella_utils)

############################################################
#                          Tests                           #
############################################################
# Compile with the following lines to run the automated tests 
#   >>  cmake . -B build -DSTELLA_ENABLE_TESTS=on
#   >>  cmake --build build -j1 --target check
# Note that the automated python tests are not working with 
# cmake. If we change the folder name back from "automated_
# stella_tests" to "integrated" then cmake does try to run 
# the tests but it will fail with the error:
#     Start 3: stella_integrated_tests
#     1/1 Test #3: stella_integrated_tests .......***Failed  
#     no such file or directory
# Since the tests are running fine with make, we will not 
# spend any time to fix this right now.
############################################################

# If -DSTELLA_ENABLE_TESTS=on we need to compile the tests as well 
if (STELLA_ENABLE_TESTS)

  # For the fortran tests, we need to download the pFUnit module
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
  
  # Include the fortran tests
  include(CTest)
  enable_testing()
  include(GKfunctions)
  set(STELLA_CTEST_CASES)
  stella_add_standard_tests(TEST_NAMES
    simple
    utils_convert
    )
  add_custom_target(run-automated-fortran-tests-with-cmake
    COMMAND ctest --output-on-failure -R unit
    DEPENDS ${STELLA_CTEST_CASES}
    )

  # With --target check we run both the fortran and python tests
  add_custom_target(check
    DEPENDS
    run-automated-fortran-tests-with-cmake
    )
else()
  set(STELLA_PFUNIT_LOCATION "Not used")
endif()

############################################################
# Configure summary

# The optimisation and/or debug flags are in the CMAKE_Fortran_FLAGS*
# variables. We need both the common flags as well as those for the
# build type actually being used. Note: this might behave weirdly on
# Windows. Might need to expand CMAKE_CONFIGURATION_TYPES instead?
############################################################

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

