# Define the wrapper library that wraps our library
add_library(${PY_NAME} SHARED bezier_com_traj.cpp)
target_link_libraries(${PY_NAME} ${PROJECT_NAME} eigenpy::eigenpy)
# don't prepend wrapper library name with lib
set_target_properties(${PY_NAME} PROPERTIES PREFIX "")

if(APPLE)
  # We need to change the extension for python bindings
  set_target_properties(${PY_NAME} PROPERTIES SUFFIX ".so")
endif(APPLE)

install(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

set(${PROJECT_NAME}_PYTHON_TESTS
    binding_tests.py
    # compare_pin_inv_dyn.py requires pinocchio_inv_dyn
)

foreach(test ${${PROJECT_NAME}_PYTHON_TESTS})
  add_python_unit_test(${test} "python/test/${test}" "python")
endforeach(test)
