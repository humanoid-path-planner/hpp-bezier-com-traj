SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
  common_solve_methods.cpp
  costfunction_definition.cpp
  solver-abstract.cpp
  eiquadprog-fast.cpp
  computeCOMTraj.cpp
  solve_0_step.cpp
  utils.cpp
  waypoints_definition.cpp
  )


if (USE_GLPK)
  SET(${LIBRARY_NAME}_SOURCES ${${LIBRARY_NAME}_SOURCES}
    glpk-wrapper.cpp
    )
endif(USE_GLPK)

ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

if (USE_GLPK)
  TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${glpk_LIBRARY})
endif(USE_GLPK)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} hpp-centroidal-dynamics)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} curves)

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
