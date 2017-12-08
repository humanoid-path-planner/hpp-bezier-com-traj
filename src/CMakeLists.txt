cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")
# no need to include directories for qpOASES as it is automatically done through pkgconfig


SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
		${INCLUDE_DIR}/bezier-com-traj/config.hh
		${INCLUDE_DIR}/bezier-com-traj/data.hh
                ${INCLUDE_DIR}/bezier-com-traj/solve.hh
                ${INCLUDE_DIR}/bezier-com-traj/common_solve_methods.hh
                ${INCLUDE_DIR}/solver/eiquadprog-fast.hpp
                solve.cpp
                common_solve_methods.cpp
                eiquadprog-fast.cpp
		#~ ${INCLUDE_DIR}/centroidal-dynamics-lib/solver_LP_abstract.hh
		#~ ${INCLUDE_DIR}/centroidal-dynamics-lib/solver_LP_qpoases.hh
		#~ ${INCLUDE_DIR}/centroidal-dynamics-lib/solver_LP_clp.hh
		#~ ${INCLUDE_DIR}/centroidal-dynamics-lib/centroidal_dynamics.hh
		#~ centroidal_dynamics.cpp
    #~ solver_LP_abstract.cpp
    #~ solver_LP_qpoases.cpp
    #~ solver_LP_clp.cpp
    #~ util.cpp
    #~ logger.cpp
    #~ stop-watch.cpp
  )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
ADD_LIBRARY(${LIBRARY_NAME} SHARED ${${LIBRARY_NAME}_SOURCES})

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

TARGET_LINK_LIBRARIES(${LIBRARY_NAME} centroidal-dynamics-lib)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} qpOASES)


SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION lib)
