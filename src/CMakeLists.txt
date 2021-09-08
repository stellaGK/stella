add_subdirectory(vmec_interface)

file(GLOB SOURCES *.f90 utils/*.f90 geo/*.f90)
file(GLOB SOURCES_FPP *.fpp utils/*.fpp )
#message(STATUS "${SOURCES}")

add_executable(stella.x  ${SOURCES} )

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message(STATUS "Found MPI Libraries in ${MPI_LIBRARIES}")
	get_filename_component (MPI_LIBRARY_PATH ${MPI_LIBRARY} DIRECTORY)
	#target_include_directories(stella.x PRIVATE ${MPI_INCLUDE_PATH} ${MPI_LIBRARY_PATH})
    target_include_directories(stella.x PRIVATE ${MPI_Fortran_INCLUDE_PATH})
	#target_link_directories(${MPI_Fortran_LIBRARIES})
	target_link_libraries( stella.x PRIVATE ${MPI_Fortran_LIBRARIES})
    #target_add_definitions(${MPI_Fortran_COMPILE_FLAGS})
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

find_package(NetCDF )

find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
	message(STATUS "Found HDF5 Libraries in ${HDF5_LIBRARIES}")
	message(STATUS "Found HDF5 Libraries in ${HDF5_INCLUDE_DIRS}")
	include_directories(${HDF5_INCLUDE_DIRS})
endif()
target_include_directories(stella.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vmec_interface/ ${CMAKE_CURRENT_SOURCE_DIR}/vmec_interface/mini_libstell)
target_include_directories(stella.x PRIVATE "/opt/netcdf-fortran-4.5.3/include")




target_link_libraries(stella.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vmec_interface/libmini_libstell.a)
#target_link_libraries(stella.x PRIVATE ${MPI_LIBRARIES})
target_link_libraries(stella.x PRIVATE "/opt/netcdf-fortran-4.5.3/lib/libnetcdff.so")

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

