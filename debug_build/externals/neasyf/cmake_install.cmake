# Install script for directory: /marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE STATIC_LIBRARY FILES "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/lib/libneasyf.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/mod/" FILES_MATCHING REGEX "/[^/]*\\.mod$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf/neasyfTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf/neasyfTargets.cmake"
         "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/CMakeFiles/Export/lib64/cmake/neasyf/neasyfTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf/neasyfTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf/neasyfTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf" TYPE FILE FILES "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/CMakeFiles/Export/lib64/cmake/neasyf/neasyfTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf" TYPE FILE FILES "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/CMakeFiles/Export/lib64/cmake/neasyf/neasyfTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/neasyf" TYPE FILE FILES
    "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/neasyfConfig.cmake"
    "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/debug_build/externals/neasyf/neasyfConfigVersion.cmake"
    "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf/cmake/GKfunctions.cmake"
    "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf/cmake/FindnetCDF.cmake"
    "/marconi_work/FUA37_STELTURB/acton/ffs_26Jul/externals/neasyf/cmake/FindnetCDFFortran.cmake"
    )
endif()

