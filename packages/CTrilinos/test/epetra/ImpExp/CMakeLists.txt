INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ImpExp_test 
  SOURCES CEpetra_ImpExp_test.cpp
  ARGS -v
  COMM mpi
  NUM_MPI_PROCS 3
  )
