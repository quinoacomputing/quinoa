
INCLUDE(TribitsAddExecutableAndTest)

# Disabled because Ifpack_CrsIct has been removed from Ifpack.
#TRIBITS_ADD_EXECUTABLE(
#  Ifpack_ex1
#  SOURCES ex1.cpp
#  ARGS
#  COMM serial mpi
#  NUM_MPI_PROCS 1-4
#  )

TRIBITS_ADD_EXECUTABLE(
  Ifpack_ex2
  SOURCES ex2.cpp
#  ARGS
  COMM serial mpi
#  NUM_MPI_PROCS 1-4
  )
