
file(GLOB_RECURSE SOURCES_LIBSTELLA *.f90 mini_libstell/*.f90 mini_libstell/*.F mini_libstell/*.f)

add_library(mini_libstell ${SOURCES_LIBSTELLA} )
target_compile_definitions(mini_libstell PRIVATE "NETCDF")
target_include_directories(mini_libstell PRIVATE "/opt/netcdf-fortran-4.5.3/include")
