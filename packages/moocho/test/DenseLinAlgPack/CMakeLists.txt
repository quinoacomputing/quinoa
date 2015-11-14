
TRIBITS_ADD_EXECUTABLE_AND_TEST( TestDenseLinAlgPack
  SOURCES
    DenseLinAlgPack_TestGenMatrixClass.cpp
    DenseLinAlgPack_TestGenMatrixOp.cpp
    DenseLinAlgPack_TestDenseLinAlgPack.cpp
    DenseLinAlgPack_TestVectorClass.cpp
    DenseLinAlgPack_TestVectorOp.cpp
    DenseLinAlgPack_TestDenseLinAlgPackIO.cpp
    TestDenseLinAlgPackMain.cpp
  NUM_MPI_PROCS 1
  )
