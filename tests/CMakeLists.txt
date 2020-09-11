SET(${PROJECT_NAME}_TESTS
  transition-quasi-static
  transition
  bezier-symbolic
  #zero-step-capturability  # hangs forever ?
  )

FOREACH(test ${${PROJECT_NAME}_TESTS})
  ADD_UNIT_TEST("${test}" "test-${test}")
  TARGET_LINK_LIBRARIES(${test} ${PROJECT_NAME})
ENDFOREACH(test)
