cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES("${INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${SPLINE_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${CDL_INCLUDE_DIR}")


SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
    ${INCLUDE_DIR}/bezier-com-traj/config.hh
    ${INCLUDE_DIR}/bezier-com-traj/data.hh
    ${INCLUDE_DIR}/bezier-com-traj/solve.hh
    ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_definition.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
    ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
    ${INCLUDE_DIR}/solver/eiquadprog-fast.hpp
    common_solve_methods.cpp
    eiquadprog-fast.cpp
    solve_transition.cpp
    solve_0_step.cpp
)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CDL_LIBRARIES})


SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
