# CMake tests specification 



TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FunctorTest
        SOURCES FunctorTest.cpp
        COMM serial mpi
        NUM_MPI_PROCS 1
)
