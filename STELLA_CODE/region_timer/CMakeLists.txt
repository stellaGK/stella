if (STELLA_ENABLE_NESMIK)
  add_library(region_timer STATIC region_timer_nesmik.c)
  target_link_libraries(region_timer PUBLIC nesmik::nesmik)
else()
  add_library(region_timer STATIC region_timer.c)
endif()

set_target_properties(region_timer PROPERTIES
  C_STANDARD 99
  C_STANDARD_REQUIRED ON
  POSITION_INDEPENDENT_CODE ON
  )

target_include_directories(region_timer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
