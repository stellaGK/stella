############################################################
#  Create the internal library <mini_libstell>
############################################################

# Source files that don't need preprocessing
set(SOURCES_LIBSTELLA_f90
  ezcdf.f90
  ezcdf_GenGet.f90
  ezcdf_inqvar.f90
  ezcdf_opncls.f90
  handle_err.f90
  makefile
  parse_extension.f
  safe_open_mod.f
  stel_constants.f
  stel_kinds.f
  vsvd0.f
  )

# Source files that need preprocessing
set(SOURCES_LIBSTELLA_fpp
  read_wout_mod.F
  )

set_source_files_properties(
  ${SOURCES_LIBSTELLA_fpp}
  PROPERTIES Fortran_PREPROCESS ON)

# Location of built and installed .mod files
set(mini_libstell_module_directory ${CMAKE_CURRENT_BINARY_DIR}/mod)
set(mini_libstell_module_install_directory "${CMAKE_INSTALL_INCLUDEDIR}" CACHE STRING
  "Fortran module installation directory")

# Create the <mini_libstell> library
add_library(mini_libstell
  ${SOURCES_LIBSTELLA_f90}
  ${SOURCES_LIBSTELLA_fpp}
  )
add_library(stella::mini_libstell ALIAS mini_libstell)
set_target_properties(mini_libstell PROPERTIES
  Fortran_MODULE_DIRECTORY "${mini_libstell_module_directory}"
  )

target_include_directories(mini_libstell PUBLIC
  $<BUILD_INTERFACE:${mini_libstell_module_directory}>
  $<INSTALL_INTERFACE:${mini_libstell_module_install_directory}>
  )

# Hard dependency on netCDF
find_package(netCDFFortran REQUIRED)
target_compile_definitions(mini_libstell PRIVATE NETCDF)
target_link_libraries(mini_libstell PUBLIC netCDF::netcdff)

# Add compiler flags for:
# - long lines
# - mismatched arguments when using the netCDF Fortran77 interface and
#   gfortran >= 10
# TODO: Fix both of these!
target_compile_options(mini_libstell PRIVATE
  $<$<Fortran_COMPILER_ID:GNU>:
    -ffree-line-length-none
    $<$<VERSION_GREATER:$<Fortran_COMPILER_VERSION>,10>:
      -fallow-argument-mismatch -w
    >
  >)

