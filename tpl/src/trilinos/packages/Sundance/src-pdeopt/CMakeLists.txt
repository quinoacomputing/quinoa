

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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Playa/src)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Core)

APPEND_SET(HEADERS
  Core/PDEOptIterCallbackBase.hpp
  Core/PDEOptLinearPDEConstrainedObj.hpp
  Core/PDEOptNonlinearPDEConstrainedObj.hpp
  Core/PDEOptPDEConstrainedObjBase.hpp
  Core/PDEOptPointData.hpp
  )

APPEND_SET(SOURCES
  Core/PDEOptIterCallbackBase.cpp
  Core/PDEOptLinearPDEConstrainedObj.cpp
  Core/PDEOptNonlinearPDEConstrainedObj.cpp
  Core/PDEOptPDEConstrainedObjBase.cpp
  Core/PDEOptPointData.cpp
  )


#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sundancePdeopt
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundanceStdFwk sundanceCore sundanceStdMesh sundanceInterop
  sundanceUtils sundancePlaya
  )
