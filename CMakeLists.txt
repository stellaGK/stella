cmake_minimum_required(VERSION 3.16)
project(stella
  VERSION 0.3
  LANGUAGES Fortran)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

#set variables here
#TODO choose which variables are necessary for stella
set(PRECISION "DOUBLE")



#set flags here
#TODO add the most common flags needed for stella

set(dialect "-ffree-line-length-none -fno-backslash -I /usr/include")
set(bounds "-fbounds-check")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} -pedantic -g -Wall -fimplicit-none")
if(PRECISION STREQUAL "DOUBLE")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()

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

add_executable(stella.x
  ${STELLA_SOURCES_f90}
  ${STELLA_SOURCES_fpp}
  )

set_target_properties(stella.x PROPERTIES
  Fortran_FORMAT FREE
)

target_compile_options(stella.x PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:
      -ffree-line-length-none>
)

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message(STATUS "Found MPI Libraries in ${MPI_LIBRARIES}")
	get_filename_component (MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)
    target_include_directories(stella.x PRIVATE ${MPI_Fortran_INCLUDE_PATH})
	target_link_libraries( stella.x PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
	message(STATUS "Found LAPACK Libraries in ${LAPACK_LIBRARIES}")
	target_link_libraries( stella.x PRIVATE ${LAPACK_LIBRARIES})
endif()


find_package(FFTW REQUIRED)
if(FFTW_FOUND)
	message(STATUS "Found FFTW Libraries in ${FFTW_LIBRARIES}")
	include_directories(${FFTW_INCLUDES})
	target_link_libraries(stella.x PRIVATE ${FFTW_LIBRARIES})
endif()


find_package(PETSc REQUIRED)
if(PETSc_FOUND)
	message(STATUS "Found PETSc Libraries in ${PETSC_LIBRARIES}")
	include_directories(${PETSC_INCLUDES})
	target_link_libraries(stella.x PRIVATE ${PETSC_LIBRARIES})
endif()


#TODO configure FindNetCDF
find_package(netCDFFortran REQUIRED)
target_link_libraries(stella.x PRIVATE netCDF::netcdff)

##################################################
# Internal libraries
add_subdirectory(geo/vmec_interface)
add_subdirectory(utils)
target_link_libraries(stella.x PRIVATE mini_libstell stella_utils)
