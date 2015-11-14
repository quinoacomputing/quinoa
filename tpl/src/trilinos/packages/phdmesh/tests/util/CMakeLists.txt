

IF(${PACKAGE_NAME}_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_util
    SOURCES 
      test_comm_all.cpp
      test_comm_bounds.cpp
      test_octtree.cpp
      test_comm.cpp
      test_containers.cpp
      test_array.cpp
      test_util.cpp
    ARGS 
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(test_util_data_copy
    SOURCE_FILES
      test_util.i
    EXEDEPS test_util
  )
ENDIF()
