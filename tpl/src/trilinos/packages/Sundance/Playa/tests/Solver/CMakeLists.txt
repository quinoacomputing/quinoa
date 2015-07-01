# CMake tests specification for Sundance/tests-solvers/SolverTests


INCLUDE(AddTestBatch)

SET(MPITests 
             PoissonTest
             BlockTriangularTest 
             PoissonBoltzmannTest )

SET(SerialOnlyTests 
             BelosPoissonTest 
             EigenTest 
	     PCGTest
             UserDefPrecondTest)


SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)

ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        )

ADD_TEST_BATCH(SerialOnlyTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )      


TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyFiles
  DEST_FILES 
                   amesos.xml 
                   anasazi-ml.xml 
                   aztec-ifpack.xml 
                   aztec-ml.xml 
                   aztec-native.xml
                   aztec.xml
                   belos-ml.xml
                   belos-ifpack.xml
                   bicgstab.xml
                   gmres.xml
                   nox.xml            
                   poissonParams.xml
                   ifpackParams.xml
                   mlParams.xml
                   userPrecParams.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialOnlyTests}
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverTestCopyMeshFiles
  DEST_FILES 
     nut-stiffness-n-2-p-1.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/Playa/etc/Meshes
  EXEDEPS ${MPITests} ${SerialOnlyTests}
  )
