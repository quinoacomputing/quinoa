function(get_mpi_compilers)

  # Find MPI
  find_package(MPI REQUIRED)

  # Find out underlying C compiler
  execute_process(
    COMMAND           ${MPI_C_COMPILER} "-showme:command"
    OUTPUT_VARIABLE   UNDERLYING_C_COMPILER
  )
  # Remove the newline character at the end of UNDERLYING_C_COMPILER
  string(REGEX REPLACE "[\r\n]" "" UNDERLYING_C_COMPILER "${UNDERLYING_C_COMPILER}")

  # Find out underlying C++ compiler
  execute_process(
    COMMAND           ${MPI_CXX_COMPILER} "-showme:command"
    OUTPUT_VARIABLE   UNDERLYING_CXX_COMPILER
  )
  # Remove the newline character at the end of UNDERLYING_CXX_COMPILER
  string(REGEX REPLACE "[\r\n]" "" UNDERLYING_CXX_COMPILER "${UNDERLYING_CXX_COMPILER}")

  # Find out underlying Fortran compiler
  execute_process(
    COMMAND           ${MPI_Fortran_COMPILER} "-showme:command"
    OUTPUT_VARIABLE   UNDERLYING_FORTRAN_COMPILER
  )
  # Remove the newline character at the end of UNDERLYING_FORTRAN_COMPILER
  string(REGEX REPLACE "[\r\n]" "" UNDERLYING_FORTRAN_COMPILER "${UNDERLYING_FORTRAN_COMPILER}")

  # Echo compilers
  MESSAGE(STATUS "MPI C compiler: " ${MPI_C_COMPILER})
  MESSAGE(STATUS "MPI C++ compiler: " ${MPI_CXX_COMPILER})
  MESSAGE(STATUS "MPI Fortran compiler: " ${MPI_Fortran_COMPILER})
  MESSAGE(STATUS "MPI underlying C compiler: " ${UNDERLYING_C_COMPILER})
  MESSAGE(STATUS "MPI underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})
  MESSAGE(STATUS "MPI underlying Fortran compiler: " ${UNDERLYING_FORTRAN_COMPILER})

endfunction()
