

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_AmesosFactory_Tridiag
  SOURCES example_AmesosFactory_Tridiag.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
