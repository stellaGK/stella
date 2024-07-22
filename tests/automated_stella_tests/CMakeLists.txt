add_test(NAME stella_integrated_tests
  COMMAND ${CMAKE_COMMAND} -E env STELLA_EXE_PATH=$<TARGET_FILE:stella> pytest -vv ${CMAKE_CURRENT_SOURCE_DIR}
  )
