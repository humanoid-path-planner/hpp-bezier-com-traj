cmake_minimum_required(VERSION 3.1)

# Project properties
set(PROJECT_NAME hpp-bezier-com-traj)
set(PROJECT_DESCRIPTION
    "Multi contact trajectory generation for the COM using Bezier curves")

# Project options
option(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
option(USE_GLPK "Use sparse lp solver" OFF)
option(PRINT_QHULL_INEQ
       "generate text file containing last inequality computed" OFF)

# Project configuration
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CUSTOM_HEADER_DIR "hpp/bezier-com-traj")
set(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
include(cmake/hpp.cmake)
include(cmake/boost.cmake)

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

if(BUILD_PYTHON_INTERFACE)
  string(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})
  add_project_dependency(eigenpy 2.7.12 REQUIRED)
endif(BUILD_PYTHON_INTERFACE)

# Project dependencies

add_project_dependency(hpp-centroidal-dynamics REQUIRED)
add_project_dependency(ndcurves 1.0.0 REQUIRED)

if(USE_GLPK)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                        "${PROJECT_SOURCE_DIR}/cmake/find-external/glpk")
  add_project_dependency(glpk REQUIRED)
  add_definitions(-DUSE_GLPK_SOLVER=1)
endif(USE_GLPK)

# Main Library
set(${PROJECT_NAME}_HEADERS
    include/${CUSTOM_HEADER_DIR}/common_solve_methods.hh
    include/${CUSTOM_HEADER_DIR}/common_solve_methods.inl
    include/${CUSTOM_HEADER_DIR}/cost/costfunction_definition.hh
    include/${CUSTOM_HEADER_DIR}/data.hh
    include/${CUSTOM_HEADER_DIR}/definitions.hh
    include/${CUSTOM_HEADER_DIR}/flags.hh
    include/${CUSTOM_HEADER_DIR}/local_config.hh
    include/${CUSTOM_HEADER_DIR}/solve_end_effector.hh
    include/${CUSTOM_HEADER_DIR}/solve.hh
    include/${CUSTOM_HEADER_DIR}/solver/eiquadprog-fast.hpp # TODO: use
                                                            # stack-of-task/eiquadprog
    include/${CUSTOM_HEADER_DIR}/solver/glpk-wrapper.hpp
    include/${CUSTOM_HEADER_DIR}/solver/solver-abstract.hpp
    include/${CUSTOM_HEADER_DIR}/utils.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_dc1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_j0_j1_ddc1_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_j0_x3_j1_ddc1_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_c0_dc0_ddc0_j0_x5_j1_ddc1_dc1_c1.hh
    include/${CUSTOM_HEADER_DIR}/waypoints/waypoints_definition.hh)

set(${PROJECT_NAME}_SOURCES
    src/common_solve_methods.cpp
    src/costfunction_definition.cpp
    src/solver-abstract.cpp
    src/eiquadprog-fast.cpp
    src/computeCOMTraj.cpp
    src/solve_0_step.cpp
    src/utils.cpp
    src/waypoints_definition.cpp)

if(USE_GLPK)
  set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES} src/glpk-wrapper.cpp)
endif(USE_GLPK)

add_library(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES}
                                   ${${PROJECT_NAME}_HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} ndcurves::ndcurves
                      hpp-centroidal-dynamics::hpp-centroidal-dynamics)

if(USE_GLPK)
  target_include_directories(${PROJECT_NAME} PUBLIC ${glpk_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} ${glpk_LIBRARY})
endif(USE_GLPK)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  DESTINATION lib)

add_subdirectory(tests)
if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)

install(FILES package.xml DESTINATION share/${PROJECT_NAME})
