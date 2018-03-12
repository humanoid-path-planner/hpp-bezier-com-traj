cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")


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

  PKG_CONFIG_USE_DEPENDENCY(${NAME} qpOASES)

  # Link against Boost and project library.
  TARGET_LINK_LIBRARIES(${NAME}
    ${Boost_LIBRARIES} centroidal-dynamics-lib bezier-com-traj
    )

ENDMACRO(ADD_TESTCASE)

ADD_TESTCASE (test-transition-quasiStatic False)



if ( MSVC )
        SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )

ADD_EXECUTABLE(bezier-com-traj-test test_static_equilibrium.cpp)

TARGET_LINK_LIBRARIES(bezier-com-traj-test bezier-com-traj)


