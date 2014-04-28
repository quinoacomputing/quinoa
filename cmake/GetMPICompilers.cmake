if(__get_mpi_compilers)
	return()
endif()
set(__get_mpi_compilers YES)


function(get_mpi_compilers)

  # Find MPI
  find_package(MPI)
  set(MPI_C_FOUND ${MPI_C_FOUND} PARENT_SCOPE)
  set(MPI_CXX_FOUND ${MPI_CXX_FOUND} PARENT_SCOPE)
  set(MPI_Fortran_FOUND ${MPI_Fortran_FOUND} PARENT_SCOPE)

  # Find out underlying C compiler
  execute_process(
    COMMAND           ${MPI_C_COMPILER} "-show"
    OUTPUT_VARIABLE   UNDERLYING_C_COMPILER
  )
  if (UNDERLYING_C_COMPILER)
    string(REPLACE " " ";" C_COMP ${UNDERLYING_C_COMPILER})
    LIST(GET C_COMP 0 UNDERLYING_C_COMPILER)
    string(REGEX REPLACE "[\r\n]" "" UNDERLYING_C_COMPILER
           "${UNDERLYING_C_COMPILER}")
    if(NOT IS_ABSOLUTE ${UNDERLYING_C_COMPILER})
      execute_process(
        COMMAND           which ${UNDERLYING_C_COMPILER}
        OUTPUT_VARIABLE   UNDERLYING_C_COMPILER
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
    set(UNDERLYING_C_COMPILER "${UNDERLYING_C_COMPILER}" PARENT_SCOPE)
  endif()

  # Find out underlying C++ compiler
  execute_process(
    COMMAND           ${MPI_CXX_COMPILER} "-show"
    OUTPUT_VARIABLE   UNDERLYING_CXX_COMPILER
  )
  if (UNDERLYING_CXX_COMPILER)
    string(REPLACE " " ";" CXX_COMP ${UNDERLYING_CXX_COMPILER})
    LIST(GET CXX_COMP 0 UNDERLYING_CXX_COMPILER)
    string(REGEX REPLACE "[\r\n]" "" UNDERLYING_CXX_COMPILER
           "${UNDERLYING_CXX_COMPILER}")
    if(NOT IS_ABSOLUTE ${UNDERLYING_CXX_COMPILER})
      execute_process(
        COMMAND           which ${UNDERLYING_CXX_COMPILER}
        OUTPUT_VARIABLE   UNDERLYING_CXX_COMPILER
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
    set(UNDERLYING_CXX_COMPILER "${UNDERLYING_CXX_COMPILER}" PARENT_SCOPE)
  endif()

  # Find out underlying Fortran compiler
  execute_process(
    COMMAND           ${MPI_Fortran_COMPILER} "-show"
    OUTPUT_VARIABLE   UNDERLYING_FORTRAN_COMPILER
  )
  if (UNDERLYING_FORTRAN_COMPILER)
    string(REPLACE " " ";" F_COMP ${UNDERLYING_FORTRAN_COMPILER})
    LIST(GET F_COMP 0 UNDERLYING_FORTRAN_COMPILER)
    string(REGEX REPLACE "[\r\n]" "" UNDERLYING_FORTRAN_COMPILER
           "${UNDERLYING_FORTRAN_COMPILER}")
    if(NOT IS_ABSOLUTE ${UNDERLYING_FORTRAN_COMPILER})
      execute_process(
        COMMAND           which ${UNDERLYING_FORTRAN_COMPILER}
        OUTPUT_VARIABLE   UNDERLYING_FORTRAN_COMPILER
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
    set(UNDERLYING_FORTRAN_COMPILER "${UNDERLYING_FORTRAN_COMPILER}" PARENT_SCOPE)
  endif()

  # Echo MPI wrappers
  if (MPI_C_COMPILER)
    MESSAGE(STATUS "MPI C wrapper: " ${MPI_C_COMPILER})
    if(NOT IS_ABSOLUTE ${MPI_C_COMPILER})
      execute_process(
        COMMAND           which ${MPI_C_COMPILER}
        OUTPUT_VARIABLE   MPI_C_COMPILER
      )
    endif()
  endif()

  if (MPI_CXX_COMPILER)
    MESSAGE(STATUS "MPI C++ wrapper: " ${MPI_CXX_COMPILER})
    if(NOT IS_ABSOLUTE ${MPI_CXX_COMPILER})
      execute_process(
        COMMAND           which ${MPI_CXX_COMPILER}
        OUTPUT_VARIABLE   MPI_CXX_COMPILER
      )
    endif()
  endif()

  if (MPI_Fortran_COMPILER)
    MESSAGE(STATUS "MPI Fortran wrapper: " ${MPI_Fortran_COMPILER})
    if(NOT IS_ABSOLUTE ${MPI_Fortran_COMPILER})
      execute_process(
        COMMAND           which ${MPI_Fortran_COMPILER}
        OUTPUT_VARIABLE   MPI_Fortran_COMPILER
      )
    endif()
  endif()

  # Echo underlying compilers
  if (UNDERLYING_C_COMPILER)
    MESSAGE(STATUS "Underlying C compiler: " ${UNDERLYING_C_COMPILER})
  endif()

  if (UNDERLYING_CXX_COMPILER)
    MESSAGE(STATUS "Underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})
  endif()

  if (UNDERLYING_FORTRAN_COMPILER)
    MESSAGE(STATUS "Underlying Fortran compiler: "
            ${UNDERLYING_FORTRAN_COMPILER})
  endif()

endfunction()
