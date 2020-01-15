cmake_minimum_required(VERSION 2.8)

SET(PROJECT_NAME hpp-bezier-com-traj)
SET(PROJECT_DESCRIPTION "Multi contact trajectory generation for the COM using Bezier curves")

INCLUDE(cmake/hpp.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)
INCLUDE(cmake/test.cmake)

SET(CUSTOM_HEADER_DIR "hpp/bezier-com-traj")

COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

# Inhibit all warning messages.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# remove flag that makes all warnings into errors
string (REPLACE "-Werror" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
OPTION (USE_GLPK "Use sparse lp solver" OFF)
OPTION (PRINT_QHULL_INEQ "generate text file containing last inequality computed" OFF)

IF(BUILD_PYTHON_INTERFACE)
  # search for python
  FINDPYTHON()
  INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS} )

  SET(BOOST_COMPONENTS python)
  SEARCH_FOR_BOOST()
  INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR} )
ENDIF(BUILD_PYTHON_INTERFACE)

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find-external/CDD")

ADD_REQUIRED_DEPENDENCY(eigen3)
ADD_REQUIRED_DEPENDENCY(hpp-centroidal-dynamics)
ADD_REQUIRED_DEPENDENCY("curves >= 0.3")
FIND_PACKAGE(CDD REQUIRED)
INCLUDE_DIRECTORIES(${CDD_INCLUDE_DIRS})

IF(USE_GLPK)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/find-external/glpk")
  find_package (glpk REQUIRED)
  add_definitions ( -DUSE_GLPK_SOLVER=1)
  INCLUDE_DIRECTORIES(${glpk_INCLUDE_DIR})
ENDIF(USE_GLPK)


add_subdirectory(include/${CUSTOM_HEADER_DIR})
add_subdirectory(src)
add_subdirectory(tests)
IF(BUILD_PYTHON_INTERFACE)
  add_subdirectory (python)
ENDIF(BUILD_PYTHON_INTERFACE)
