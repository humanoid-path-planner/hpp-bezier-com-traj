cmake_minimum_required( VERSION 2.8 )

include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")

ADD_REQUIRED_DEPENDENCY("eigenpy")


# Define the wrapper library that wraps our library
add_library( bezier_com_traj SHARED bezier_com_traj.cpp )
target_link_libraries( bezier_com_traj ${Boost_LIBRARIES} ${PROJECT_NAME} )
# don't prepend wrapper library name with lib
set_target_properties( bezier_com_traj PROPERTIES PREFIX "" )

IF(APPLE)
	# We need to change the extension for python bindings
	SET_TARGET_PROPERTIES(bezier_com_traj PROPERTIES SUFFIX ".so")
ENDIF(APPLE)
PKG_CONFIG_USE_DEPENDENCY(bezier_com_traj eigenpy)

INSTALL(
	FILES ${LIBRARY_OUTPUT_PATH}/bezier_com_traj.so	DESTINATION ${PYTHON_SITELIB}
)
