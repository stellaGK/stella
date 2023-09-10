
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was neasyfConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../../../../../usr/local" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# If using the build directory directly, we need the CMake modules too
if(EXISTS "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf/cmake")
  list(APPEND CMAKE_MODULE_PATH "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf/cmake")
endif()

include(GKfunctions)

if(EXISTS "")
  set(netCDFFortran_ROOT "")
endif()

find_dependency(netCDFFortran 4.4.4)

include("${CMAKE_CURRENT_LIST_DIR}/neasyfTargets.cmake")
