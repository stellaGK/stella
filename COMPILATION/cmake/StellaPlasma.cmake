############################################################
#  PLASMA linear algebra library (STELLA_ENABLE_PLASMA)
############################################################
# Provides the variable STELLA_PLASMA_LIBRARIES, which lists the
# targets that need to be linked to use PLASMA.
#
# An installed PLASMA is used if one is found, in this order:
#  1. find_package(plasma CONFIG), i.e. the plasmaConfig.cmake
#     that PLASMA installs, found through PLASMA_DIR (cmake or
#     environment variable, the install prefix of PLASMA),
#     plasma_DIR, plasma_ROOT or CMAKE_PREFIX_PATH;
#  2. pkg-config (plasma.pc, which PLASMA also installs);
#  3. the libraries themselves, in PLASMA_DIR or the default
#     library paths.
# Otherwise the git submodule EXTERNALS/plasma is compiled as an
# external project in the build directory. Set
# STELLA_PLASMA_FROM_SUBMODULE=ON to always use the submodule.
#
# PLASMA requires the CBLAS and LAPACKE libraries, extra cmake
# flags for the submodule can be passed through PLASMA_CMAKE_FLAGS
# (a cmake list). Only the PLASMA libraries are compiled, since
# the PLASMA tester needs LAPACKE routines missing from older
# LAPACK versions, e.g. LAPACKE_zlangb in LAPACK 3.10 (Ubuntu
# 22.04). The libraries are used from the PLASMA build directory,
# since installing PLASMA would also install (and therefore
# require) the tester.
############################################################

option(STELLA_PLASMA_FROM_SUBMODULE "Compile PLASMA from the EXTERNALS/plasma submodule, even if PLASMA is installed" OFF)

if (NOT STELLA_ENABLE_DOUBLE)
  message(FATAL_ERROR "STELLA_ENABLE_PLASMA requires STELLA_ENABLE_DOUBLE, since stella is bound to the double precision PLASMA routines")
endif()

if (NOT PLASMA_DIR AND DEFINED ENV{PLASMA_DIR})
  set(PLASMA_DIR "$ENV{PLASMA_DIR}")
endif()

set(STELLA_PLASMA_LIBRARIES)

if (NOT STELLA_PLASMA_FROM_SUBMODULE)

  # 1. The CMake package of an installed PLASMA. Its plasmaConfig.cmake does
  # not include the file that defines the targets, so we include it ourselves
  set(_plasma_hints)
  if (PLASMA_DIR)
    list(APPEND _plasma_hints "${PLASMA_DIR}" "${PLASMA_DIR}/lib/cmake/plasma" "${PLASMA_DIR}/lib64/cmake/plasma")
  endif()
  find_package(plasma CONFIG QUIET HINTS ${_plasma_hints})
  if (plasma_FOUND AND NOT TARGET plasma::plasma AND EXISTS "${plasma_DIR}/plasmaTargets.cmake")
    include("${plasma_DIR}/plasmaTargets.cmake")
  endif()
  if (TARGET plasma::plasma AND TARGET plasma::plasma_core_blas)
    message(STATUS "Using the installed PLASMA from ${plasma_DIR}")
    set(STELLA_PLASMA_LIBRARIES plasma::plasma plasma::plasma_core_blas)
  endif()

  # 2. pkg-config
  if (NOT STELLA_PLASMA_LIBRARIES)
    find_package(PkgConfig QUIET)
    if (PKG_CONFIG_FOUND)
      pkg_check_modules(STELLA_PKG_PLASMA QUIET IMPORTED_TARGET plasma)
      if (STELLA_PKG_PLASMA_FOUND)
        message(STATUS "Using the installed PLASMA from pkg-config (${STELLA_PKG_PLASMA_LIBDIR})")
        set(STELLA_PLASMA_LIBRARIES PkgConfig::STELLA_PKG_PLASMA)
      endif()
    endif()
  endif()

  # 3. The libraries themselves
  if (NOT STELLA_PLASMA_LIBRARIES)
    if (PLASMA_DIR)
      set(_plasma_find_options PATHS "${PLASMA_DIR}" PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
    endif()
    find_library(STELLA_PLASMA_LIBRARY plasma ${_plasma_find_options})
    find_library(STELLA_PLASMA_CORE_BLAS_LIBRARY plasma_core_blas ${_plasma_find_options})
    if (STELLA_PLASMA_LIBRARY AND STELLA_PLASMA_CORE_BLAS_LIBRARY)
      message(STATUS "Using the installed PLASMA libraries ${STELLA_PLASMA_LIBRARY}")
      set(STELLA_PLASMA_LIBRARIES ${STELLA_PLASMA_LIBRARY} ${STELLA_PLASMA_CORE_BLAS_LIBRARY})
    elseif (PLASMA_DIR)
      message(FATAL_ERROR "Could not find the PLASMA libraries in PLASMA_DIR=${PLASMA_DIR}")
    endif()
  endif()

endif()

# 4. Compile the git submodule
if (NOT STELLA_PLASMA_LIBRARIES)

  get_filename_component(_plasma_source_dir "${CMAKE_CURRENT_LIST_DIR}/../../EXTERNALS/plasma" ABSOLUTE)
  set(_plasma_build_dir "${CMAKE_CURRENT_BINARY_DIR}/EXTERNALS/plasma/build")
  if (NOT EXISTS "${_plasma_source_dir}/CMakeLists.txt")
    message(FATAL_ERROR "No installed PLASMA was found and the PLASMA submodule is missing, run: "
      "git submodule update --init EXTERNALS/plasma, or build without PLASMA: -DSTELLA_ENABLE_PLASMA=OFF")
  endif()
  message(STATUS "No installed PLASMA found, compiling PLASMA from the submodule ${_plasma_source_dir}")

  set(_plasma_library "${_plasma_build_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}plasma${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(_plasma_core_blas_library "${_plasma_build_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}plasma_core_blas${CMAKE_SHARED_LIBRARY_SUFFIX}")

  include(ExternalProject)
  ExternalProject_Add(plasma_external
    SOURCE_DIR "${_plasma_source_dir}"
    BINARY_DIR "${_plasma_build_dir}"
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Release
      ${PLASMA_CMAKE_FLAGS}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target plasma plasma_core_blas
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS "${_plasma_library}" "${_plasma_core_blas_library}"
    )

  # Imported targets for the PLASMA libraries
  add_library(stella_plasma_core_blas SHARED IMPORTED GLOBAL)
  set_target_properties(stella_plasma_core_blas PROPERTIES IMPORTED_LOCATION "${_plasma_core_blas_library}")
  add_library(stella_plasma SHARED IMPORTED GLOBAL)
  set_target_properties(stella_plasma PROPERTIES IMPORTED_LOCATION "${_plasma_library}")
  add_dependencies(stella_plasma_core_blas plasma_external)
  add_dependencies(stella_plasma plasma_external)
  set(STELLA_PLASMA_LIBRARIES stella_plasma stella_plasma_core_blas)

endif()
