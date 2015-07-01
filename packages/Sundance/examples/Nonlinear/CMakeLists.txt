
INCLUDE(AddTestBatch)


SET(MPITests
)

SET(SerialTests
  FixedPointBratu1D
  HandLinearizedNewtonBratu1D
  AutoLinearizedNewtonBratu1D
  FullyAutomatedNewtonBratu1D
  ContinuationBratu1D
)

ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        NUM_MPI_PROCS 4
        )


ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(NonlinearExampleCopyFiles
  DEST_FILES 
                   amesos.xml 
                   playa-newton-amesos.xml 
                   nox-amesos.xml 
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/examples/Nonlinear
  EXEDEPS ${MPITests} ${SerialTests}
  )



