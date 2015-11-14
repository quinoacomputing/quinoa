
ADD_SUBDIRECTORIES(
  RowMatrix
  ImpExp
  )


TRIBITS_INCLUDE_DIRECTORIES( REQUIRED_DURING_INSTALLATION_TESTING
                     ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_SOURCE_DIR}/.. )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Comm_UnitTests
  SOURCES
    CEpetra_Comm_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_SerialComm_UnitTests
  SOURCES
    CEpetra_SerialComm_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_MpiComm_UnitTests
  SOURCES
    CEpetra_MpiComm_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_BlockMap_UnitTests
  SOURCES
    CEpetra_BlockMap_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Map_UnitTests
  SOURCES
    CEpetra_Map_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Object_UnitTests
  SOURCES
    CEpetra_Object_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Distributor_UnitTests
  SOURCES
    CEpetra_Distributor_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Directory_UnitTests
  SOURCES
    CEpetra_Directory_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Import_UnitTests
  SOURCES
    CEpetra_Import_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Export_UnitTests
  SOURCES
    CEpetra_Export_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  CEpetra_OffsetIndex_UnitTests
#  SOURCES
#    CEpetra_OffsetIndex_UnitTests.cpp
#    ../CTrilinos_UnitTestHelpers.cpp
#    ${TEUCHOS_STD_UNIT_TEST_MAIN}
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_DistObject_UnitTests
  SOURCES
    CEpetra_DistObject_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_SrcDistObject_UnitTests
  SOURCES
    CEpetra_SrcDistObject_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_BLAS_UnitTests
  SOURCES
    CEpetra_BLAS_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Flops_UnitTests
  SOURCES
    CEpetra_Flops_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_CompObject_UnitTests
  SOURCES
    CEpetra_CompObject_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_MultiVector_UnitTests
  SOURCES
    CEpetra_MultiVector_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Vector_UnitTests
  SOURCES
    CEpetra_Vector_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_CrsGraph_UnitTests
  SOURCES
    CEpetra_CrsGraph_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_CrsMatrix_UnitTests
  SOURCES
    CEpetra_CrsMatrix_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Operator_UnitTests
  SOURCES
    CEpetra_Operator_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_RowMatrix_UnitTests
  SOURCES
    CEpetra_RowMatrix_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_JadMatrix_UnitTests
  SOURCES
    CEpetra_JadMatrix_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_Time_UnitTests
  SOURCES
    CEpetra_Time_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_LinearProblem_UnitTests
  SOURCES
    CEpetra_LinearProblem_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CEpetra_LAPACK_UnitTests
  SOURCES
    CEpetra_LAPACK_UnitTests.cpp
    ../CTrilinos_UnitTestHelpers.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

