# - Returns the underlying C compiler for mpicc

if(__get_underlying_compiler)
	return()
endif()
set(__get_underlying_compiler YES)

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

# Echo compilers
MESSAGE(STATUS "MPI C compiler: " ${MPI_C_COMPILER})
MESSAGE(STATUS "MPI C++ compiler: " ${MPI_CXX_COMPILER})
MESSAGE(STATUS "MPI underlying C compiler: " ${UNDERLYING_C_COMPILER})
MESSAGE(STATUS "MPI underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})
