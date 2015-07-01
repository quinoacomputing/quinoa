# CMake tests specification for Sundance/tests-solvers/OperatorTests


INCLUDE(AddTestBatch)

SET(SerialTests
  UnconstrainedTests)


ADD_TEST_BATCH(
        SerialTests
        COMM serial mpi
        NUM_MPI_PROCS 1
)

SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/Playa/etc/SolverParameters)

TRIBITS_COPY_FILES_TO_BINARY_DIR(OptimizerTestCopyFiles
  DEST_FILES 
  steepestDescent.xml 
  basicLMBFGS.xml 
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )
