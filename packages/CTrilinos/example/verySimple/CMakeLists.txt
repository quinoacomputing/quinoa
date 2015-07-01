
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  verySimple
  SOURCES verySimple.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

