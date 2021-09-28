# Add an alias for an imported target
# Workaround for CMake < 3.18
# Taken from https://github.com/conan-io/conan/issues/2125#issuecomment-351176653
function(gk_add_library_alias dst src)
  add_library(${dst} INTERFACE IMPORTED)
  foreach(name INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS INTERFACE_COMPILE_OPTIONS)
    get_property(value TARGET ${src} PROPERTY ${name} )
    set_property(TARGET ${dst} PROPERTY ${name} ${value})
  endforeach()
endfunction()


# Call nx-config with an argument, and append the resulting path to a list
# Taken from https://github.com/LiamBindle/geos-chem/blob/feature/CMake/CMakeScripts/FindNetCDF.cmake
function(gk_inspect_netcdf_config VAR NX_CONFIG ARG)
  execute_process(
    COMMAND ${NX_CONFIG} ${ARG}
    OUTPUT_VARIABLE NX_CONFIG_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NX_CONFIG_OUTPUT)
    set(${VAR} ${NX_CONFIG_OUTPUT} PARENT_SCOPE)
  endif()
endfunction()

# Adds target for pfunit tests in test_source
# name utils_tests_${test_name} and adds to
# list of know tests cases to give to ctest.
function(gk_add_test test_source test_name)
  add_pfunit_test(utils_tests_${test_name}
    "${test_source}" "" "")
  target_link_libraries(utils_tests_${test_name} gkutils)
  list(APPEND GK_CTEST_CASES utils_tests_${test_name})
  set(GK_CTEST_CASES ${GK_CTEST_CASES} PARENT_SCOPE)
endfunction()

# Helper function to easily add multiple separate
# tests provided they exist at tests/test_${name}.pf
# and we're happy to identify them as utils_tests_${name}
function(gk_add_standard_tests)
  cmake_parse_arguments(
        GK_ADD "" "" "TEST_NAMES" ${ARGN}
    )
  foreach(name ${GK_ADD_TEST_NAMES})
    gk_add_test("tests/test_${name}.pf" ${name})
  endforeach()
  set(GK_CTEST_CASES ${GK_CTEST_CASES} PARENT_SCOPE)
endfunction()
