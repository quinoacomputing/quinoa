
INCLUDE(AddTestBatch)

SET(MPITests
  Skeleton
  CommandLineOptions
)

SET(SerialTests
  Laplace3D
  XMLParameterList
)




ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(GettingStartedCopyFiles
  DEST_FILES 
     plateWithHole3D-1.exo
     paramExample.xml
     aztec-ml.xml
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/examples/GettingStarted
  EXEDEPS ${MPITests} ${SerialTests}
  )



