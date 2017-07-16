################################################################################
#
# \file      cmake/MPICompilers.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
# \brief     Find the MPI wrappers
#
################################################################################

if(__get_mpi_compilers)
	return()
endif()
set(__get_mpi_compilers YES)

macro(get_compiler CMAKE_MPI_COMPILER UNDERLYING_COMPILER)

  set(UNDERLYING_COMPILER "")

  # Detect a Cray machine. This is based on testing for the following
  # environment variables. At least one of these is always defined on a Cray,
  # depending on what programming environment is loaded among the modules.
  if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
      NOT DEFINED ENV{CRAY_PRGENVGNU} AND
      NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
      NOT DEFINED ENV{CRAY_PRGENVINTEL})

    execute_process(
      COMMAND           ${CMAKE_MPI_COMPILER} "-show"
      OUTPUT_VARIABLE   TMP_COMPILER
      ERROR_VARIABLE    ERR
    )
    if(ERR)
      MESSAGE(FATAL_ERROR
              "Executing command '${CMAKE_MPI_COMPILER} -show' failed.")
    endif()

    if (TMP_COMPILER)
      string(REPLACE " " ";" COMP ${TMP_COMPILER})
      LIST(GET COMP 0 TMP_COMPILER)
      string(REGEX REPLACE "[\r\n]" "" TMP_COMPILER "${TMP_COMPILER}")
      if(NOT IS_ABSOLUTE ${TMP_COMPILER})
        execute_process(
          COMMAND           which ${TMP_COMPILER}
          OUTPUT_VARIABLE   TMP_COMPILER
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      endif()
      set(${UNDERLYING_COMPILER} "${TMP_COMPILER}")
      set(${UNDERLYING_COMPILER} "${TMP_COMPILER}" PARENT_SCOPE)
    endif()
  else()
    # For Cray machines, the underlying compiler is the same as the wrapper;
    # export it to both local and parent scope.
    message(STATUS "Cray detected, not querying underlying compiler for "
                   ${CMAKE_MPI_COMPILER})
    set(${UNDERLYING_COMPILER} "${CMAKE_MPI_COMPILER}")
    set(${UNDERLYING_COMPILER} "${CMAKE_MPI_COMPILER}" PARENT_SCOPE)
  endif()

endmacro()


function(get_mpi_compilers)

  # Find MPI
  find_package(MPI REQUIRED)
  # Export to parent scope
  set(MPI_C_FOUND ${MPI_C_FOUND} PARENT_SCOPE)
  set(MPI_CXX_FOUND ${MPI_CXX_FOUND} PARENT_SCOPE)
  set(MPI_Fortran_FOUND ${MPI_Fortran_FOUND} PARENT_SCOPE)

  # Find underlying C, C++, and Fortran compilers
  if (MPI_C_COMPILER)
    get_compiler(${MPI_C_COMPILER} UNDERLYING_C_COMPILER)
  endif()

  if (MPI_CXX_COMPILER)
    get_compiler(${MPI_CXX_COMPILER} UNDERLYING_CXX_COMPILER)
  endif()

  if (MPI_Fortran_COMPILER)
    get_compiler(${MPI_Fortran_COMPILER} UNDERLYING_Fortran_COMPILER)
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

  # Echo underling compilers
  if (UNDERLYING_C_COMPILER)
    MESSAGE(STATUS "Underlying C compiler: " ${UNDERLYING_C_COMPILER})
  endif()

  if (UNDERLYING_CXX_COMPILER)
    MESSAGE(STATUS "Underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})
  endif()

  if (UNDERLYING_Fortran_COMPILER)
    MESSAGE(STATUS "Underlying Fortran compiler: " ${UNDERLYING_Fortran_COMPILER})
  endif()

endfunction()
