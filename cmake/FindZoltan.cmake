# Find the Zoltan library from Sandia
#
#  ZOLTAN_FOUND - System has both ExodusII and Nemesis
#  ZOLTAN_INCLUDES - The ExodusII include directory
#  ZOLTAN_LIBRARIES - The libraries needed to use ExodusII
#
#  Set ZOLTAN_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(ZOLTAN_ROOT "/path/to/custom/zoltan") # prefer over system
#  find_package(Zoltan)
#  if(ZOLTAN_FOUND)
#    target_link_libraries (TARGET ${ZOLTAN_LIBRARIES})
#  endif()

if(ZOLTAN_INCLUDES AND ZOLTAN_LIBRARIES)
  # Already in cache, be silent
  set (ZOLTAN_FIND_QUIETLY TRUE)
endif()

FIND_PATH(ZOLTAN_INCLUDES NAMES zoltan.h HINTS ${ZOLTAN_ROOT}/include)
FIND_LIBRARY(ZOLTAN_LIBRARIES NAMES zoltan HINTS ${ZOLTAN_ROOT}/lib)

# Handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Zoltan DEFAULT_MSG ZOLTAN_LIBRARIES ZOLTAN_INCLUDES)

MARK_AS_ADVANCED(ZOLTAN_INCLUDES ZOLTAN_LIBRARIES)
