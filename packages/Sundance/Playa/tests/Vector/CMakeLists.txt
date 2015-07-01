# CMake tests specification for Playa/tests/Vectors


INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
        BasicVectorTest
        SOURCES BasicVectorTest.cpp
        COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
        BlockVectorTest
        SOURCES BlockVectorTest.cpp
        COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LCTest
  SOURCES LCTest.cpp
  COMM serial mpi
)
