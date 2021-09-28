cmake_minimum_required( VERSION 3.5)
project(stella)
enable_language(Fortran)



set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

#set variables here
#TODO choose which variables are necessary for stella
set(PRECISION "DOUBLE")

add_subdirectory(geo/vmec_interface)

add_executable(stella.x  ${SOURCES} )

#set flags here
#TODO add the most common flags needed for stella

set(dialect "-ffree-line-length-none -fno-backslash -I /usr/include")
set(bounds "-fbounds-check")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} -pedantic -g -Wall -fimplicit-none")
if(PRECISION STREQUAL "DOUBLE")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()

#get package includes and libs here

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

#link executable to internal library
target_include_directories(stella.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vmec_interface/)
target_link_libraries(stella.x PRIVATE mini_libstell)


#TODO configure FindNetCDF
find_package(netCDFFortran REQUIRED)
target_link_libraries(stella.x PRIVATE netCDF::netcdff)


#get all fpp and f90 files
file(GLOB SOURCES *.f90 utils/*.f90 geo/*.f90)
file(GLOB SOURCES_FPP *.fpp utils/*.fpp )

#generate f90 files from fpp files and link them to executable
set( GENERATED_SOURCES "" )
foreach ( SOURCE_FPP ${SOURCES_FPP} ) 
    get_filename_component (SOURCE_NO_EXT ${SOURCE_FPP} NAME_WE )
    set(SOURCE_F90 "${SOURCE_NO_EXT}.f90" )
	message(${SOURCE_F90})
	add_custom_command(
    	OUTPUT  "${SOURCE_F90}"
    	COMMAND cpp -P -traditional -DFCOMPILER=_GFORTRAN_ -ffreestanding  -DISO_C_BINDING -DMPI -DFFT=_FFTW3_ -DNETCDF -DLAPACK ${SOURCE_FPP} ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_F90}
    	DEPENDS ${SOURCE_FPP}
	)
    target_sources( stella.x PRIVATE ${SOURCE_F90} )
endforeach()
