# CMake tests specification 


INCLUDE(AddTestBatch)


SET(Tests 
        Circle2D
)
    
ADD_TEST_BATCH(Tests 
        COMM serial mpi
        NUM_MPI_PROCS 1
        )



