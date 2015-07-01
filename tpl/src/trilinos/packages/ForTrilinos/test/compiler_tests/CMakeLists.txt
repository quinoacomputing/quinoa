INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  interoperability_check
  SOURCES
    check_compiler_for_C_interoperability.F90
    ${PACKAGE_SOURCE_DIR}/test/compiler_tests/check_compiler_for_C_interoperability.F90
  COMM serial mpi
  NUM_MPI_PROCS 1
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

IF(${PACKAGE_NAME}_ENABLE_OBJECT_ORIENTED)
 
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
  check_for_memory_leaks
  SOURCES 
    check_compiler_for_memory_leaks.F90
    ${PACKAGE_SOURCE_DIR}/test/compiler_tests/check_compiler_for_memory_leaks.F90
  COMM serial mpi
  NUM_MPI_PROCS 4
  ARGS -v
  STANDARD_PASS_OUTPUT
 )

 #TRIBITS_ADD_EXECUTABLE_AND_TEST(
 #  check_compiler_for_leaks_with_arrays
 # SOURCES 
 #   check_compiler_for_leaks_with_arrays.F90
 #   ${PACKAGE_SOURCE_DIR}/test/compiler_tests/check_compiler_for_leaks_with_arrays.F90
 # COMM serial mpi
 # NUM_MPI_PROCS 4
 # ARGS -v
 # STANDARD_PASS_OUTPUT
 #)

ENDIF()
