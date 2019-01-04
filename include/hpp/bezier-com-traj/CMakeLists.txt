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
  cost/costfunction_definition.hh
  waypoints/waypoints_definition.hh
  waypoints/waypoints_c0_dc0_c1.hh
  waypoints/waypoints_c0_dc0_dc1.hh
  waypoints/waypoints_c0_dc0_dc1_c1.hh
  waypoints/waypoints_c0_dc0_ddc0_c1.hh
  waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
  waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
  include/solver/solver-abstract.hpp
  include/solver/eiquadprog-fast.hpp

  )

IF(USE_GLPK)
  SET(${PROJECT_NAME}_HEADERS ${${PROJECT_NAME}_HEADERS}
    include/solver/glpk-wrapper.hpp
    )
ENDIF(USE_GLPK)
