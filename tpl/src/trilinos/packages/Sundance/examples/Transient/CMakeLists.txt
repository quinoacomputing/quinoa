
INCLUDE(AddTestBatch)

SET(SerialTests
  TransientHeat2D
  StepControlledTransientNonlinear1D
)

SET(MPITests
  TransientNonlinear1D
)




ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TransientExampleCopyFiles
  DEST_FILES 
                   amesos.xml 
                   aztec-ml.xml 
                   playa-newton-amesos.xml
                   playa-newton-aztec-ml.xml
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/examples/Transient
  EXEDEPS ${MPITests} ${SerialTests}
  )



