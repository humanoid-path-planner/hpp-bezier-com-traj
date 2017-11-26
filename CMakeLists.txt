cmake_minimum_required(VERSION 2.6)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/python.cmake)
INCLUDE(cmake2/FindqpOASES.cmake)

SET(PROJECT_NAME bezier-com-traj)
SET(PROJECT_DESCRIPTION
	 "Multi contact trajectory generation for the COM using Bezier curves"
 )
SET(PROJECT_URL "")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${PROJECT_SOURCE_DIR}/cmake2/")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(EXTERNAL_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/external/lib")

if ( MSVC )  
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")
endif()

SETUP_PROJECT()

# Inhibit all warning messages.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# remove flag that makes all warnings into errors
string (REPLACE "-Werror" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" OFF)
IF(BUILD_PYTHON_INTERFACE)
# search for python
	FINDPYTHON(2.7  REQUIRED)
	find_package( PythonLibs 2.7 REQUIRED )
	include_directories( ${PYTHON_INCLUDE_DIRS} )

	find_package( Boost COMPONENTS python REQUIRED )
	include_directories( ${Boost_INCLUDE_DIR} )

ENDIF(BUILD_PYTHON_INTERFACE)


# Declare Headers
SET(${PROJECT_NAME}_HEADERS
                include/bezier-com-traj/config.hh
                include/bezier-com-traj/solve.hh
                include/bezier-com-traj/common_solve_methods.hh
                include/solver/eiquadprog-fast.hpp
  )

find_package(Eigen3)
if (EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3")
else()
        message(STATUS "Eigen3 not found looking for Eigen2")
        find_package(Eigen2 REQUIRED)
endif()


#SEARCH_FOR_QPOASES()
ADD_REQUIRED_DEPENDENCY("qpOASES")

add_subdirectory (src)
add_subdirectory (test)
#IF(BUILD_PYTHON_INTERFACE)
#	add_subdirectory (python)
#ENDIF(BUILD_PYTHON_INTERFACE)

SETUP_PROJECT_FINALIZE()

