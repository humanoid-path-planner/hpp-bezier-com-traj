cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES("${INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${SPLINE_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${CDL_INCLUDE_DIR}")
if (USE_GLPK)
INCLUDE_DIRECTORIES("${GLPK_INCLUDE_DIR}")
endif(USE_GLPK)


SET(LIBRARY_NAME ${PROJECT_NAME})

if (USE_GLPK)
    SET(${LIBRARY_NAME}_SOURCES
        ${INCLUDE_DIR}/bezier-com-traj/config.hh
        ${INCLUDE_DIR}/bezier-com-traj/definitions.hh
        ${INCLUDE_DIR}/bezier-com-traj/data.hh
        ${INCLUDE_DIR}/bezier-com-traj/utils.hh
        ${INCLUDE_DIR}/bezier-com-traj/flags.hh
        ${INCLUDE_DIR}/bezier-com-traj/solve.hh
        ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.hh
        ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.inl
        ${INCLUDE_DIR}/bezier-com-traj/cost/costfunction_definition.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_definition.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
        ${INCLUDE_DIR}/solver/eiquadprog-fast.hpp
        ${INCLUDE_DIR}/solver/glpk-wrapper.hpp
        ${INCLUDE_DIR}/solver/solver-abstract.hpp
        common_solve_methods.cpp
        costfunction_definition.cpp
        solver-abstract.cpp
        eiquadprog-fast.cpp
        glpk-wrapper.cpp
        computeCOMTraj.cpp
        solve_0_step.cpp
        utils.cpp
        waypoints_definition.cpp
    )
else()
    SET(${LIBRARY_NAME}_SOURCES
        ${INCLUDE_DIR}/bezier-com-traj/config.hh
        ${INCLUDE_DIR}/bezier-com-traj/definitions.hh
        ${INCLUDE_DIR}/bezier-com-traj/data.hh
        ${INCLUDE_DIR}/bezier-com-traj/utils.hh
        ${INCLUDE_DIR}/bezier-com-traj/flags.hh
        ${INCLUDE_DIR}/bezier-com-traj/solve.hh
        ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.hh
        ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.inl
        ${INCLUDE_DIR}/bezier-com-traj/cost/costfunction_definition.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_definition.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
        ${INCLUDE_DIR}/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
        ${INCLUDE_DIR}/solver/eiquadprog-fast.hpp
        ${INCLUDE_DIR}/solver/solver-abstract.hpp
        common_solve_methods.cpp
        costfunction_definition.cpp
        solver-abstract.cpp
        eiquadprog-fast.cpp
        computeCOMTraj.cpp
        solve_0_step.cpp
        utils.cpp
        waypoints_definition.cpp
    )
endif(USE_GLPK)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

if (USE_GLPK)
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${GLPK_LIB})
endif(USE_GLPK)
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CDL_LIBRARIES})


SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
