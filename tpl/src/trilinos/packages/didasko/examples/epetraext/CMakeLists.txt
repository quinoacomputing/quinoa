
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EpetraEXT_ex3
  SOURCES ex3.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  )
