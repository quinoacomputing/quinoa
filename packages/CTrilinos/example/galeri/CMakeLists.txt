
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  example_Galeri_LinearProblem
  SOURCES example_Galeri_LinearProblem.c
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
