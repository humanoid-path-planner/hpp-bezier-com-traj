cmake_minimum_required(VERSION 2.6)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/python.cmake)

SET(PROJECT_NAME bezier-com-traj)
SET(PROJECT_DESCRIPTION
	 "Multi contact trajectory generation for the COM using Bezier curves"
 )
SET(PROJECT_URL "")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${PROJECT_SOURCE_DIR}/cmake2/")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(EXTERNAL_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/external/lib")

if ( MSVC )  
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")
endif()


# Where to compile shared objects
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)


SETUP_PROJECT()

# Inhibit all warning messages.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# remove flag that makes all warnings into errors
string (REPLACE "-Werror" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
OPTION (USE_GLPK "Use sparse lp solver" ON)

IF(BUILD_PYTHON_INTERFACE)
# search for python
	FINDPYTHON(2.7  REQUIRED)
	find_package( PythonLibs 2.7 REQUIRED )
	include_directories( ${PYTHON_INCLUDE_DIRS} )

	find_package( Boost COMPONENTS python REQUIRED )
	include_directories( ${Boost_INCLUDE_DIR} )

ENDIF(BUILD_PYTHON_INTERFACE)


find_package (centroidal-dynamics-lib REQUIRED)
IF(USE_GLPK)
find_package (glpk REQUIRED)
add_definitions ( -DUSE_GLPK_SOLVER=1)
ENDIF(USE_GLPK)
find_package (Spline REQUIRED)

# Declare Headers
IF(USE_GLPK)
SET(${PROJECT_NAME}_HEADERS
                include/bezier-com-traj/data.hh
                include/bezier-com-traj/utils.hh
                include/bezier-com-traj/flags.hh
                include/bezier-com-traj/definitions.hh
                include/bezier-com-traj/config.hh
                include/bezier-com-traj/solve.hh
                include/bezier-com-traj/solve_end_effector.hh
                include/bezier-com-traj/common_solve_methods.hh
                include/bezier-com-traj/common_solve_methods.inl
                include/bezier-com-traj/cost/costfunction_definition.hh
                include/bezier-com-traj/waypoints/waypoints_definition.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
                include/solver/solver-abstract.hpp
                include/solver/glpk-wrapper.hpp
                include/solver/eiquadprog-fast.hpp

  )
else()
SET(${PROJECT_NAME}_HEADERS
                include/bezier-com-traj/data.hh
                include/bezier-com-traj/utils.hh
                include/bezier-com-traj/flags.hh
                include/bezier-com-traj/definitions.hh
                include/bezier-com-traj/config.hh
                include/bezier-com-traj/solve.hh
                include/bezier-com-traj/solve_end_effector.hh
                include/bezier-com-traj/common_solve_methods.hh
                include/bezier-com-traj/common_solve_methods.inl
                include/bezier-com-traj/cost/costfunction_definition.hh
                include/bezier-com-traj/waypoints/waypoints_definition.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_dc1_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_dc1_c1.hh
                include/bezier-com-traj/waypoints/waypoints_c0_dc0_ddc0_ddc1_dc1_c1.hh
                include/solver/solver-abstract.hpp
                include/solver/eiquadprog-fast.hpp

  )
ENDIF(USE_GLPK)

find_package(Eigen3)
if (EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3")
else()
        message(STATUS "Eigen3 not found looking for Eigen2")
        find_package(Eigen2 REQUIRED)
endif()


add_subdirectory(src)
add_subdirectory(tests)
IF(BUILD_PYTHON_INTERFACE)
        add_subdirectory (python)
ENDIF(BUILD_PYTHON_INTERFACE)

SETUP_PROJECT_FINALIZE()

