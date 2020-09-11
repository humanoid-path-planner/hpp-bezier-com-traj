cmake_minimum_required(VERSION 3.1)

# Project properties
SET(PROJECT_NAME hpp-bezier-com-traj)
SET(PROJECT_DESCRIPTION "Multi contact trajectory generation for the COM using Bezier curves")

# Project options
OPTION(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
OPTION(USE_GLPK "Use sparse lp solver" OFF)
OPTION(PRINT_QHULL_INEQ "generate text file containing last inequality computed" OFF)

# Project configuration
SET(PROJECT_USE_CMAKE_EXPORT TRUE)
SET(CUSTOM_HEADER_DIR "hpp/bezier-com-traj")
SET(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
INCLUDE(cmake/hpp.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)

# Project definition
COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

IF(BUILD_PYTHON_INTERFACE)
  FINDPYTHON()
  FIND_NUMPY()
  STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  SEARCH_FOR_BOOST_PYTHON(REQUIRED)
  ADD_PROJECT_DEPENDENCY(eigenpy REQUIRED)
ENDIF(BUILD_PYTHON_INTERFACE)

# Project dependencies

ADD_PROJECT_DEPENDENCY(hpp-centroidal-dynamics REQUIRED)
ADD_PROJECT_DEPENDENCY(curves REQUIRED)

IF(USE_GLPK)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/find-external/glpk")
  ADD_PROJECT_DEPENDENCY(glpk REQUIRED)
  add_definitions ( -DUSE_GLPK_SOLVER=1)
ENDIF(USE_GLPK)

# Main Library
SET(${PROJECT_NAME}_HEADERS
  include/${CUSTOM_HEADER_DIR}/data.hh
  include/${CUSTOM_HEADER_DIR}/utils.hh
  include/${CUSTOM_HEADER_DIR}/flags.hh
  include/${CUSTOM_HEADER_DIR}/definitions.hh
  include/${CUSTOM_HEADER_DIR}/local_config.hh
  include/${CUSTOM_HEADER_DIR}/solve.hh
  include/${CUSTOM_HEADER_DIR}/solve_end_effector.hh
  include/${CUSTOM_HEADER_DIR}/common_solve_methods.hh
  include/${CUSTOM_HEADER_DIR}/common_solve_methods.inl
  include/${CUSTOM_HEADER_DIR}/cost/costfunction_definition.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_definition.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_c1.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_dc1.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_dc1_c1.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_c1.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
  include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
  include/${CUSTOM_HEADER_DIR}/solver/solver-abstract.hpp
  include/${CUSTOM_HEADER_DIR}/solver/eiquadprog-fast.hpp  # TODO: use stack-of-task/eiquadprog
  )

SET(${PROJECT_NAME}_SOURCES
  src/common_solve_methods.cpp
  src/costfunction_definition.cpp
  src/solver-abstract.cpp
  src/eiquadprog-fast.cpp
  src/computeCOMTraj.cpp
  src/solve_0_step.cpp
  src/utils.cpp
  src/waypoints_definition.cpp
  )

if(USE_GLPK)
  SET(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES}
    glpk-wrapper.cpp
    )
endif(USE_GLPK)

ADD_LIBRARY(${PROJECT_NAME} SHARED
  ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} curves::curves hpp-centroidal-dynamics::hpp-centroidal-dynamics)

if(USE_GLPK)
  TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} ${glpk_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${glpk_LIBRARY})
endif(USE_GLPK)

INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)

ADD_SUBDIRECTORY(tests)
IF(BUILD_PYTHON_INTERFACE)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)
