SET(${PROJECT_NAME}_WAYPOINTS_HEADERS
  waypoints_definition.hh
  waypoints_c0_dc0_c1.hh
  waypoints_c0_dc0_dc1.hh
  waypoints_c0_dc0_dc1_c1.hh
  waypoints_c0_dc0_ddc0_c1.hh
  waypoints_c0_dc0_ddc0_dc1_c1.hh
  waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
  )

INSTALL(FILES
  ${${PROJECT_NAME}_WAYPOINTS_HEADERS}
  DESTINATION include/${CUSTOM_HEADER_DIR}/waypoints
  )
