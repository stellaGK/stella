cmake_minimum_required( VERSION 3.5)
project(stella)
enable_language(Fortran)



set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

#set variables here
#TODO choose which variables are necessary for stella
set(PRECISION "DOUBLE")




add_subdirectory(${PROJECT_SOURCE_DIR}/src)
