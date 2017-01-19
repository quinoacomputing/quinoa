################################################################################
#
# \file      cmake/FindZoltan2.cmake
# \author    C. Junghans
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Zoltan2 header
# \date      Thu 19 Jan 2017 10:32:46 AM MST
#
################################################################################

# Find Zoltan2 headers and libraries
#
#  ZOLTAN2_FOUND        - True if Zoltan2 is found
#  ZOLTAN2_INCLUDE_DIRS - Zoltan2 include files paths
#  ZOLTAN2_LIBRARIES    - List of Zoltan2 libraries
#
#  Set ZOLTAN2_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(ZOLTAN2_ROOT "/path/to/custom/zoltan2") # prefer over system
#  find_package(ZOLTAN2)
#  if(ZOLTAN2_FOUND)
#    target_link_libraries (TARGET ${ZOLTAN2_LIBRARIES})
#  endif()

# If already in cache, be silent
if (ZOLTAN2_INCLUDE_DIRS AND ZOLTAN2_LIBRARIES)
  set (ZOLTAN2_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(ZOLTAN2_INCLUDE_DIR NAMES Zoltan2_PartitioningSolution.hpp
          PATH_SUFFIXES trilinos
          HINTS ${ZOLTAN2_ROOT}/include $ENV{ZOLTAN2_ROOT}/include)

# Look for the library
if(NOT BUILD_SHARED_LIBS)
  FIND_LIBRARY(ZOLTAN2_LIBRARY NAMES libzoltan2.a libtrilinos_zoltan2.a
               PATH_SUFFIXES trilinos
               HINTS ${ZOLTAN2_ROOT}/lib $ENV{ZOLTAN2_ROOT}/lib)
else()
  FIND_LIBRARY(ZOLTAN2_LIBRARY NAMES zoltan2 trilinos_zoltan2
               PATH_SUFFIXES trilinos
               HINTS ${ZOLTAN2_ROOT}/lib $ENV{ZOLTAN2_ROOT}/lib)
endif()

set(ZOLTAN2_LIBRARIES ${ZOLTAN2_LIBRARY} )
set(ZOLTAN2_INCLUDE_DIRS ${ZOLTAN2_INCLUDE_DIR} )

# Handle the QUIETLY and REQUIRED arguments and set ZOLTAN2_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Zoltan2 DEFAULT_MSG ZOLTAN2_INCLUDE_DIRS ZOLTAN2_LIBRARIES)

MARK_AS_ADVANCED(ZOLTAN2_INCLUDE_DIRS ZOLTAN2_LIBRARIES)
