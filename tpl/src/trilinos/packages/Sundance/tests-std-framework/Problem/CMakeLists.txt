

INCLUDE(AddTestBatch)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

SET(DisabledTests
  FunctionalTest3D 
  transientOneZoneEddy2D
)

SET(SerialTests
  BigCellFilterTest
  BesselTest2D 
  PCDNavierStokesCouette2D
  NonlinearPDECOTest1D
  NonlinearPDECOTest2D
  FunctionalTest 
  PerforatedPlate
  PerforatedPlate_QuadInt
  Stokes2D_HN_Q1
  PoissonPeano2D
  PoissonPeano3D
  StokesPeano2D 
  MixedTest
  RefinementTest2D
  Poisson2D_HN_Nitsch
  Stokes_Chanel_Nitsche_line_obstic
  TransientNonlinTest
  ControlledTransient1D
  TriBdryTest
)



SET(MPITests 
             LPTests1D
             Poisson2D
             CrankNicolsonHeat1D
             FieldIOTest
             AdvectionDiffusion2D 
             LinearVariationalTest 
             OptPoisson1D
             VorticityNavierStokes2D
             SteadyBurgers1D
             MixedTest_HN
             NavStok_Chanel_stat_HN
             Poisson3DTest_HN_Q1Q1
             Poisson3DTest_HN_Q2Q1
             LinearPartialDomain_mixed
             NonlinearPartialDomain_mixed
             NavStok_Chanel_stat_HN_Nitsch
             Poisson3D_Surf
   )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MPIProbTestBatch
    SOURCES ProbMPITestDriver.cpp 
       BlockStochPoissonTest1D.cpp 
       HighOrderPoisson2D.cpp 
       HighOrderPoissonBernstein2D.cpp 
       HighOrderProjection2D.cpp 
       NonlinReducedIntegration.cpp
       SpectralPoisson1D.cpp 
       SpectralSqrt.cpp
    COMM serial mpi
    NUM_MPI_PROCS 4
  )



TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SerialProbTestBatch
    SOURCES ProbSerialTestDriver.cpp 
       NonlinearPartialDomain.cpp
       NonlinearPeriodic1D.cpp
       LinearPeriodic1D.cpp
       PoissonOnDisk.cpp
       TetQuadTransformationTest.cpp
       DiscFunc3D.cpp
       Kepler.cpp
       CNBugTest.cpp
       EdgeDFTest.cpp
       AToCDensitySample.cpp
       DuffingFloquet.cpp
       SecondOrderFloquet.cpp
       SubmaximalDFTest.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
  )


SET(SolverParamPath ${PACKAGE_SOURCE_DIR}/etc/SolverParameters)


ADD_TEST_BATCH(MPITests 
        COMM serial mpi
        )

ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )

TRIBITS_COPY_FILES_TO_BINARY_DIR(TestCopyFiles
  DEST_FILES 
                   amesos.xml 
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
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS ${MPITests} ${SerialTests} ${SerialProbsTestBatch} ${MPIProbsTestBatch}
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(TestMeshCopyFiles
  DEST_FILES 
    couette-2.exo
    couette-3.exo
    disk.ncdf
    disk.1.ele
    disk.1.node
    diskWithHole2D-0.exo
    diskWithHole2D-1.exo
    tetPrism.exo
    square.1.ele
    square.1.edge
    square.1.node
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/tests-std-framework/Problem
  EXEDEPS ${MPITests} ${SerialTests} ${SerialProbsTestBatch} ${MPIProbsTestBatch}
  )



# Eikonal2D needs special handling because of the fortran file

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Eikonal2D
       SOURCES Eikonal2D.cpp bessi0.f
       COMM serial mpi
       NUM_MPI_PROCS 1
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(EikonalCopyFiles
  DEST_FILES nox.xml
  SOURCE_DIR ${SolverParamPath}
  EXEDEPS Eikonal2D
  )
