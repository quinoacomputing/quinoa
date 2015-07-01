
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  Triutils_ex1
  SOURCES ex1.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )

TRIBITS_ADD_EXECUTABLE(
  Triutils_ex2
  SOURCES ex2.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )
