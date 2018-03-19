cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/tests)
INCLUDE_DIRECTORIES("${INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${SPLINE_INCLUDE_DIR}")
INCLUDE_DIRECTORIES("${CDL_INCLUDE_DIR}")


# Make Boost.Test generates the main function in test cases.
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)
# ADD_TESTCASE(NAME)
# ------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(ADD_TESTCASE NAME GENERATED)
  IF (${GENERATED} STREQUAL TRUE)
    ADD_EXECUTABLE(${NAME} ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cc)
  ELSE()
    ADD_EXECUTABLE(${NAME} ${NAME}.cc)
  ENDIF()
  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})

  # Link against Boost and project library.
  TARGET_LINK_LIBRARIES(${NAME}
    ${Boost_LIBRARIES} bezier-com-traj  ${CDL_LIBRARIES}
    )

ENDMACRO(ADD_TESTCASE)

ADD_TESTCASE (test-transition-quasiStatic False)
ADD_TESTCASE (test-transition False)

if ( MSVC )
        SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

ADD_EXECUTABLE(bezier-com-traj-test test-zero-step-capturability.cpp test_helper.hh)

TARGET_LINK_LIBRARIES(bezier-com-traj-test bezier-com-traj ${CDL_LIBRARIES})


