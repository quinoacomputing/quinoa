
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex1
  SOURCES ex1.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex2
  SOURCES ex2.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex3
  SOURCES ex3.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex4
  SOURCES ex4.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex5
  SOURCES ex5.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex6
  SOURCES ex6.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Teuchos_ex7
  SOURCES ex7.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )
