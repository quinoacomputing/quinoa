

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(SundanceConfig.h)
TRIBITS_CONFIGURE_FILE(SundanceDefaultPath.hpp)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/SundanceDefaultPath.hpp
  )

#
# Utils files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Utilities)

APPEND_SET(HEADERS
  Utilities/SundanceADReal.hpp
  Utilities/SundanceArrayOfTuples.hpp
  Utilities/SundanceCollectiveExceptionCheck.hpp
  Utilities/SundanceCombinatorialUtils.hpp
  Utilities/SundanceDebug.hpp
  Utilities/SundanceDefs.hpp
  Utilities/SundanceEnumTypeField.hpp
  Utilities/SundanceFunctorDomain.hpp
  Utilities/SundanceIntHashSet.hpp
  Utilities/SundanceIntVec.hpp
  Utilities/SundanceMap.hpp
  Utilities/SundanceMultiSet.hpp
  Utilities/SundanceNamedObject.hpp
  Utilities/SundanceNoncopyable.hpp
  Utilities/SundanceObjectWithInstanceID.hpp
  Utilities/SundanceObjectWithVerbosity.hpp
  Utilities/SundanceOrderedHandle.hpp
  Utilities/SundanceOrderedTuple.hpp
  Utilities/SundanceOut.hpp
  Utilities/SundanceParamUtils.hpp
  Utilities/SundancePathUtils.hpp
  Utilities/SundancePoint.hpp
  Utilities/SundanceSet.hpp
  Utilities/SundanceStdMathFunctors.hpp
  Utilities/SundanceTypeUtils.hpp
  Utilities/SundanceUnaryFunctor.hpp
  Utilities/SundanceVersionString.hpp
  Utilities/SundanceVertexSort.hpp
  Utilities/SundanceWatchFlag.hpp
  )

APPEND_SET(SOURCES
  Utilities/SundanceADReal.cpp
  Utilities/SundanceCollectiveExceptionCheck.cpp
  Utilities/SundanceCombinatorialUtils.cpp
  Utilities/SundanceFunctorDomain.cpp
  Utilities/SundanceIntHashSet.cpp
  Utilities/SundanceIntVec.cpp
  Utilities/SundanceParamUtils.cpp
  Utilities/SundancePathUtils.cpp
  Utilities/SundancePoint.cpp
  Utilities/SundanceStdMathFunctors.cpp
  Utilities/SundanceUnaryFunctor.cpp
  Utilities/SundanceVertexSort.cpp
  Utilities/SundanceWatchFlag.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  sundanceUtils
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS sundancePlaya
  )
