

TRIBITS_INCLUDE_DIRECTORIES( REQUIRED_DURING_INSTALLATION_TESTING
                     ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CAmesos_UnitTests
  SOURCES
    CAmesos_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  CAmesos_BaseSolver_UnitTests
#  SOURCES
#    CAmesos_BaseSolver_UnitTests.cpp
#    ../CTrilinos_UnitTestHelpers.cpp
#    ${TEUCHOS_STD_UNIT_TEST_MAIN}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )

