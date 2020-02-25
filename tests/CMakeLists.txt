SET(${PROJECT_NAME}_TESTS
  transition-quasi-static
  transition
  bezier-symbolic
  #zero-step-capturability  # hangs forever ?
  )

FOREACH(test ${${PROJECT_NAME}_TESTS})
  ADD_UNIT_TEST("${test}" "test-${test}")
  PKG_CONFIG_USE_DEPENDENCY(${test} hpp-centroidal-dynamics)
  PKG_CONFIG_USE_DEPENDENCY(${test} curves)
  PKG_CONFIG_USE_DEPENDENCY(${test} eigen3)
  TARGET_LINK_LIBRARIES(${test} ${Boost_LIBRARIES} ${PROJECT_NAME} hpp-centroidal-dynamics)
ENDFOREACH(test)
