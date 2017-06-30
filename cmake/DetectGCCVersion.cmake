################################################################################
#
# \file      cmake/DetectGCCVersion.cmake
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Detect GNU C version
#
################################################################################

execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                OUTPUT_VARIABLE GCC_VERSION)
string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})

list(LENGTH GCC_VERSION_COMPONENTS GCCVER)
if (GCCVER GREATER 0)
  list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
endif()
if (GCCVER GREATER 1)
  list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
endif()

message(STATUS "GCC version: ${GCC_MAJOR}.${GCC_MINOR}")
