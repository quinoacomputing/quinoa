

INCLUDE(AddTestBatch)


SET(SerialTests
  NonlinControl1D
  AlgebraicConstrainedOpt
)



SET(PDEOptNLSolverParamPath ${PACKAGE_SOURCE_DIR}/Playa/etc/SolverParameters)


ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(PDEOptNLTestCopyFiles
  DEST_FILES 
                   amesos.xml 
                   basicLMBFGS.xml
                   anasazi-ml.xml 
                   aztec-ifpack.xml 
                   aztec-ml.xml 
                   aztec-native.xml
                   aztec.xml
                   belos-ifpack.xml 
                   belos-ml.xml 
                   bicgstab.xml
                   gmres.xml
                   nox.xml            
                   nox-amesos.xml            
                   nox-aztec.xml      
                   nox-newton.xml   
                   pcd-outer-belos.xml               
                   poissonParams.xml
                   userPrecParams.xml
  SOURCE_DIR ${PDEOptNLSolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests}
  )


