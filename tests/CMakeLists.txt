set(${PROJECT_NAME}_TESTS transition-quasi-static transition bezier-symbolic
                          # zero-step-capturability  # hangs forever ?
)

foreach(test ${${PROJECT_NAME}_TESTS})
  add_unit_test("${test}" "test-${test}.cpp")
  target_link_libraries(${test} ${PROJECT_NAME})
endforeach(test)
