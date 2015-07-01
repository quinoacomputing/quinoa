
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ADT_3D_Burgers_6th_Pade 
  SOURCES main.F90 vector_field.F90 scalar_field.F90 initializer.F90 field.F90 globals.F90
  COMM serial mpi
  NUM_MPI_PROCS 4
  ARGS -v
  STANDARD_PASS_OUTPUT
  )

#SET_SOURCE_FILES_PROPERTIES( ${SOURCES} PROPERTIES LANGUAGE Fortran )

