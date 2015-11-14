

#
# A) Package-specific configuration options
#


#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#
# Stubs files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-utils/Utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Stubs)

APPEND_SET(HEADERS
  Stubs/SundanceCellFilterStub.hpp
  Stubs/SundanceNullCellFilterStub.hpp
  )

APPEND_SET(SOURCES
  Stubs/SundanceCellFilterStub.cpp
  Stubs/SundanceNullCellFilterStub.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sundanceInterop
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundanceUtils
  )
