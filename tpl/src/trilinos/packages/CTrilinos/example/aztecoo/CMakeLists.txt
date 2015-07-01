
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  simple_epetra_aztecoo
  SOURCES simple_epetra_aztecoo.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
