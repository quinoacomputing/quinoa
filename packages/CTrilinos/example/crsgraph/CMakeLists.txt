
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsGraph_example
  SOURCES CrsGraph_example.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  XHOSTTYPE Windows Windows_NT
  )

