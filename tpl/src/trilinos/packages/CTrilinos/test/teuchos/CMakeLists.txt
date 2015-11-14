

TRIBITS_INCLUDE_DIRECTORIES( REQUIRED_DURING_INSTALLATION_TESTING
                     ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_any_UnitTests
  SOURCES
    CTeuchos_any_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_CommandLineProcessor_UnitTests
  SOURCES
    CTeuchos_CommandLineProcessor_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_ParameterEntry_UnitTests
  SOURCES
    CTeuchos_ParameterEntry_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTeuchos_ParameterList_UnitTests
  SOURCES
    CTeuchos_ParameterList_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

