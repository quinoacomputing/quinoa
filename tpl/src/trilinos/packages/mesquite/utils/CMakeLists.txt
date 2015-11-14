
IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")
ENDIF()

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

APPEND_SET(HEADERS
  Mesquite_CLArgs.hpp
  Mesquite_CLArgFlag.hpp
  Mesquite_ManPage.hpp
  Mesquite_domain.hpp
  )

APPEND_SET(SOURCES
  Mesquite_CLArgs.cpp
  Mesquite_CLArgFlag.cpp
  Mesquite_ManPage.cpp
  Mesquite_domain.cpp
  )

IF (TPL_ENABLE_MPI)

ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (NOT ${PROJECT_NAME} STREQUAL "Mesquite")

#  ADD_LIBRARY( msqutil STATIC ${SOURCES} ${HEADERS} )
  TRIBITS_ADD_LIBRARY(
    msqutil
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS mesquite
    )

  TRIBITS_ADD_EXECUTABLE(
    vtkrandom
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES Mesquite_randomize.cpp
    DEPLIBS msqutil
    )

  TRIBITS_ADD_EXECUTABLE(
    msqshape
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES Mesquite_shape-improve.cpp
    DEPLIBS msqutil
    )

  TRIBITS_ADD_EXECUTABLE(
    vtkxform
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES Mesquite_transform.cpp
    DEPLIBS msqutil
    )

  TRIBITS_ADD_EXECUTABLE(
    msqquality
    NOEXEPREFIX
    NOEXESUFFIX
    SOURCES Mesquite_quality.cpp Mesquite_domain.cpp
    DEPLIBS msqutil
    )

ELSE()

  ADD_LIBRARY( msqutil STATIC ${SOURCES} ${HEADERS} )
  ADD_EXECUTABLE( vtkrandom Mesquite_randomize.cpp Mesquite_domain.cpp )
  ADD_EXECUTABLE( msqshape  Mesquite_shape-improve.cpp Mesquite_domain.cpp )
  ADD_EXECUTABLE( vtkxform  Mesquite_transform.cpp )
  ADD_EXECUTABLE( msqquality Mesquite_quality.cpp )
  TARGET_LINK_LIBRARIES( vtkrandom  msqutil mesquite )
  TARGET_LINK_LIBRARIES( msqshape   msqutil mesquite )
  TARGET_LINK_LIBRARIES( vtkxform   msqutil mesquite )
  TARGET_LINK_LIBRARIES( msqquality msqutil mesquite )

ENDIF()

IF(${PACKAGE_NAME}_ENABLE_SHARED)
  SET( MESQUITE_UTIL_DEFINES "${MESQUITE_DEFINES}" )
ELSE()
  SET( MESQUITE_UTIL_DEFINES "${MESQUITE_DEFINES} -DMESQUITE_STATIC_LIB" )
ENDIF()
SET(UTIL_SOURCES
  Mesquite_CLArgFlag.cpp 
  Mesquite_CLArgs.cpp 
  Mesquite_domain.cpp 
  Mesquite_ManPage.cpp 
  Mesquite_quality.cpp 
  Mesquite_randomize.cpp
  Mesquite_shape-improve.cpp
  Mesquite_transform.cpp
  )
SET_SOURCE_FILES_PROPERTIES( ${UTIL_SOURCES} COMPILE_FLAGS ${MESQUITE_UTIL_DEFINES} )
