# Define the wrapper library that wraps our library
add_library(${PY_NAME} SHARED bezier_com_traj.cpp )
target_link_libraries(${PY_NAME} ${PROJECT_NAME} eigenpy::eigenpy)
# don't prepend wrapper library name with lib
set_target_properties(${PY_NAME} PROPERTIES PREFIX "")

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)

INSTALL(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

SET(${PROJECT_NAME}_PYTHON_TESTS
  binding_tests.py
  #compare_pin_inv_dyn.py requires pinocchio_inv_dyn
  )

FOREACH(test ${${PROJECT_NAME}_PYTHON_TESTS})
  ADD_PYTHON_UNIT_TEST(${test} "python/test/${test}" "python")
ENDFOREACH(test)
