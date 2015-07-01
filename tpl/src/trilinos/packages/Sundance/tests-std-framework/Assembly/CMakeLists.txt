# CMake tests specification for Sundance/tests-solvers/SolverTests


INCLUDE(AddTestBatch)

SET(SerialTests BasisCheck TransformedIntegral2D  QuadratureTest RTDOFTest)


ADD_TEST_BATCH(SerialTests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )      



