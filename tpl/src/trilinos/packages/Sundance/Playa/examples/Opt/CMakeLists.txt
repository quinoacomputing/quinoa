# 


INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PlayaOptRosenbrock
  SOURCES PlayaOptRosenbrock.cpp PlayaRosenbrock.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/Playa/etc/SolverParameters)

TRIBITS_COPY_FILES_TO_BINARY_DIR(OptimizerExampleCopyFiles
  DEST_FILES 
  basicLMBFGS.xml 
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )

