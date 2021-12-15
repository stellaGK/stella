cmake_minimum_required(VERSION 3.18)

project(fortran-git
  DESCRIPTION "Interface to get git version/commit in compiled code"
  LANGUAGES Fortran
  )

add_library(fortran_git
  src/git_version.f90
  src/git_version_impl.fpp)

add_library(fortran_git::fortran_git ALIAS fortran_git)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(GetGitRevisionDescription)
git_get_commit(GIT_SHA1)
message(STATUS "GIT_SHA1=${GIT_SHA1}")
git_get_version(GIT_VERSION)
message(STATUS "GIT_VERSION=${GIT_VERSION}")
git_local_changes(GIT_STATE)
message(STATUS "GIT_STATE=${GIT_STATE}")
git_get_commit_date(GIT_DATE)
message(STATUS "GIT_DATE=${GIT_DATE}")

set_source_files_properties(git_version_impl.F90
  PROPERTIES Fortran_PREPROCESS ON)

target_compile_definitions(fortran_git PRIVATE
  "GIT_SHA1='${GIT_SHA1}'"
  "GIT_VERSION='${GIT_VERSION}'"
  "GIT_STATE='${GIT_STATE}'"
  "GIT_DATE='${GIT_DATE}'"
  )

set(fortran_git_module_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(fortran_git PROPERTIES
  Fortran_MODULE_DIRECTORY ${fortran_git_module_dir}
  Fortran_FORMAT FREE
  )

target_include_directories(fortran_git PUBLIC
  $<BUILD_INTERFACE:${fortran_git_module_dir}>
  )

##################################################
# Example use

string(COMPARE EQUAL
  "${PROJECT_NAME}" "${CMAKE_PROJECT_NAME}"
  PROJECT_IS_TOP_LEVEL
)

option(FORTRAN_GIT_BUILD_EXAMPLES "Build the example" ${PROJECT_IS_TOP_LEVEL})

if (FORTRAN_GIT_BUILD_EXAMPLES)
  add_executable(example example/example.f90)
  target_link_libraries(example fortran_git)
endif()
