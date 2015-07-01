INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST( TestIterationPack
  SOURCES
    IterationPack_AlgorithmStepTesting.cpp
    IterationPack_AlgorithmTrackTesting.cpp
    IterationPack_TestAlgorithm.cpp
    IterationPack_TestAlgorithmState.cpp
    IterationPack_TestIterationPack.cpp
    IterationPack_TestIterQuantityAccessContinuous.cpp
    TestIterationPackMain.cpp
  NUM_MPI_PROCS 1
  )
