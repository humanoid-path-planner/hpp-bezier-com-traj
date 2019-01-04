SET(${PROJECT_NAME}_HEADERS
  data.hh
  utils.hh
  flags.hh
  definitions.hh
  local_config.hh
  solve.hh
  solve_end_effector.hh
  common_solve_methods.hh
  common_solve_methods.inl
  )

INSTALL(FILES
  ${${PROJECT_NAME}_HEADERS}
  DESTINATION include/${CUSTOM_HEADER_DIR}
  )

ADD_SUBDIRECTORY(cost)
ADD_SUBDIRECTORY(waypoints)
ADD_SUBDIRECTORY(solver)
