
INCLUDE(AddTestBatch)

SET(SerialTests
  CoolingIdealFlow
  Coupled2D
  NitschePoisson2D 
  PCDStokesCouette2D
)

SET(MPITests
  ParameterSweep
  TaylorHoodStokes2D
  VorticityStokes2D 
)



ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        NUM_MPI_PROCS 4
        )


ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(LinearExampleCopyFiles
  DEST_FILES 
                   amesos.xml 
                   aztec-ifpack.xml 
                   aztec-ml.xml 
                   aztec-native.xml
                   belos-ml.xml 
                   bicgstab.xml 
                   pcd-outer-belos.xml 
                   pipe2D-1.exo
                   pipe2D-2.exo
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/examples/Linear
  EXEDEPS ${MPITests} ${SerialTests}
  )



