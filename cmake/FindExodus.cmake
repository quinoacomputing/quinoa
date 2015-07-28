# Find the Exodus II and Nemesis finite element data model library from Sandia
#
#  EXODUS_FOUND - System has both ExodusII and Nemesis libraries and exodiff
#  EXODUS_INCLUDES - The ExodusII include directory
#  EXODUS_LIBRARIES - The libraries needed to use ExodusII
#  NEMESIS_INCLUDES - The Nemesis include directory
#  NEMESIS_LIBRARIES - The libraries needed to use Nemesis
#  EXODIFF_EXECUTABLE - Exodiff used to diff two ExodusII databases
#
#  Set EXODUS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(EXODUS_ROOT "/path/to/custom/exodus_and_nemesis") # prefer over system
#  find_package(Exodus)
#  if(EXODUS_FOUND)
#    target_link_libraries (TARGET ${EXODUS_LIBRARIES})
#    target_link_libraries (TARGET ${NEMESIS_LIBRARIES})
#  endif()

if(EXODUS_INCLUDES AND EXODUS_LIBRARIES AND NEMESIS_INCLUDES AND NEMESISLIBRARIES AND EXODIFF_EXECUTABLE)
  # Already in cache, be silent
  set (EXODUS_FIND_QUIETLY TRUE)
endif()

FIND_PATH(EXODUS_INCLUDES NAMES exodusII.h HINTS ${EXODUS_ROOT}/include)
FIND_LIBRARY(EXODUS_LIBRARIES NAMES exodus exoIIv2 HINTS ${EXODUS_ROOT}/lib)

FIND_PATH(NEMESIS_INCLUDES NAMES ne_nemesisI.h HINTS ${EXODUS_ROOT}/include)
FIND_LIBRARY(NEMESIS_LIBRARIES NAMES nemesis HINTS ${EXODUS_ROOT}/lib)

FIND_PATH(EXODIFF_EXECUTABLE NAMES exodiff HINTS ${EXODUS_ROOT}/bin)

# handle the QUIETLY and REQUIRED arguments and set EXODUS_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Exodus DEFAULT_MSG EXODUS_LIBRARIES EXODUS_INCLUDES NEMESIS_LIBRARIES NEMESIS_INCLUDES EXODIFF_EXECUTABLE)

MARK_AS_ADVANCED(EXODUS_INCLUDES EXODUS_LIBRARIES NEMESIS_INCLUDES NEMESIS_LIBRARIES EXODIFF_EXECUTABLE)
