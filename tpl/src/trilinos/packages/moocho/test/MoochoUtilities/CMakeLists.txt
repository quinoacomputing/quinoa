INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST( test_stop_watch
  SOURCES TestStopWatchMain.cpp
  NUM_MPI_PROCS 1
  )
