cmake_minimum_required( VERSION 3.5)
project(stella)
enable_language(Fortran)


set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
#
##file(STRINGS myfile.txt myfile)
#message(STATUS "Content in ${myfile}")
#set(OPT "False")
#set(USE_MPI "False")
#set(USE_FFT "False")
#set(USE_NETCDF "False")
#set(USE_HDF5 "False")
#set(USE_LAPACK "False")
set(PRECISION "DOUBLE")
#
##include(myfile.txt)
#


set(dialect "-ffree-line-length-none -fno-backslash -I /usr/include")
set(bounds "-fbounds-check")
#
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} -march=native")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} -pedantic -g -Wall -fimplicit-none")
#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS}")
if(PRECISION STREQUAL "DOUBLE")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
endif()




#if(OPT)
#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
#endif()
#if(USE_OPENMP)
#	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
#endif()
#
#
#

#
#find_package(LAPACK REQUIRED)
#if(LAPACK_FOUND)
#	message(STATUS "Found LAPACK Libraries in ${LAPACK_LIBRARIES}")
#	include_directories(${LAPACK_LIBRARIES})
#endif()

#set (NETCDF_F90 "YES")
#find_package(PETSc REQUIRED)
#find_package(FFTW REQUIRED)


#TODO
#FIND NETCDF
#FIND PETSc

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
