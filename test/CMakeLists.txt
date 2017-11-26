cmake_minimum_required(VERSION 2.6)

include_directories("${SRC_DIR}")
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")

PROJECT(bezier-com-traj-test)

if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

add_executable(bezier-com-traj-test test_static_equilibrium.cpp)

TARGET_LINK_LIBRARIES(bezier-com-traj-test bezier-com-traj)
#~ TARGET_LINK_LIBRARIES(polytopetest polytope ${SRC_DIR}/../external/cddlib-094b/lib-src/libcdd.a)
