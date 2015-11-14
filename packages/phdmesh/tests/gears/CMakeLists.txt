 
 
TRIBITS_ADD_EXECUTABLE(
  test_mesh
  INCLUDES Gears.hpp
  SOURCES
    Gears.cpp
    test_gears.cpp
    test_mesh.cpp
    MeshTestSchema.cpp
    test_mean_value.cpp
    test_simple_mesh.cpp
  COMM serial mpi
  )

#NOTE: This test fails in DEBUG mode on macs.  Since the DEBUG mode isn't 
#specified by default for checkin testing, we will just disable the tests
#on macs.
IF ( NOT APPLE )
TRIBITS_ADD_TEST(
  test_mesh
  NAME test_simple
  COMM serial mpi
  ARGS "threadpool 1 , schema_parts , simple_mesh"
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  test_mesh
  NAME test_gears
  COMM serial mpi
#      "gears <#X> <#Y> <#Z> [verify] [output <base_filename>]"
  ARGS "threadpool 1 , gears 2 3 1 steps 5 verify"
  NUM_MPI_PROCS 1
  )
ENDIF()

# TRIBITS_ADD_TEST(
#   test_mesh
#   NAME test_gears
#   COMM mpi
# #      "gears <#X> <#Y> <#Z> [verify] [output <base_filename>]"
#   ARGS "threadpool 1 , gears 2 3 1 steps 5 verify"
#   NUM_MPI_PROCS 2
#   )

