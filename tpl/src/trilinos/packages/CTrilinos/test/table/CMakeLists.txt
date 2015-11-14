

TRIBITS_INCLUDE_DIRECTORIES( REQUIRED_DURING_INSTALLATION_TESTING
                     ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_Table_UnitTests
  SOURCES
    CTrilinos_Table_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_TableRepos_UnitTests
  SOURCES
    CTrilinos_TableRepos_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CTrilinos_utils_UnitTests
  SOURCES
    CTrilinos_utils_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  CTrilinos_Table_memory
#  SOURCES
#    CTrilinos_Table_memory.c
#  ARGS 12
#  COMM serial
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )

