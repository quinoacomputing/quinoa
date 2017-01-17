################################################################################
#
# \file      cmake/FindZoltan2.cmake
# \author    C. Junghans
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Zoltan2 header
# \date      Tue 17 Jan 2017 09:52:44 AM MST
#
################################################################################

# Find Zoltan2 headers and libraries
#
#  Zoltan2_FOUND        - True if Zoltan2 is found
#  Zoltan2_INCLUDE_DIRS - Zoltan2 include files paths
#  Zoltan2_LIBRARIES    - List of Zoltan2 libraries
#
#  Set Zoltan2_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(Zoltan2_ROOT "/path/to/custom/boost/mpl") # prefer over system
#  find_package(Zoltan2)
#  if(Zoltan2_FOUND)
#    target_link_libraries (TARGET ${Zoltan2_LIBRARIES})
#  endif()

# If already in cache, be silent
if (Zoltan2_INCLUDE_DIR AND Zoltan2_LIBRARY)
  set (Zoltan2_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(Zoltan2_INCLUDE_DIR NAMES Zoltan2_PartitioningSolution.hpp
          PATH_SUFFIXES trilinos
          HINTS ${Zoltan2_ROOT}/include)

# Look for the library
if(NOT BUILD_SHARED_LIBS)
  FIND_LIBRARY(Zoltan2_LIBRARY NAMES libzoltan2.a libtrilinos_zoltan2.a
               PATH_SUFFIXES trilinos
               HINTS ${Zoltan2_ROOT}/lib)
else()
  FIND_LIBRARY(Zoltan2_LIBRARY NAMES zoltan2 trilinos_zoltan2
               PATH_SUFFIXES trilinos
               HINTS ${Zoltan2_ROOT}/lib)
endif()

set(Zoltan2_LIBRARIES ${Zoltan2_LIBRARY} )
set(Zoltan2_INCLUDE_DIRS ${Zoltan2_INCLUDE_DIR} )

# Handle the QUIETLY and REQUIRED arguments and set Zoltan2_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Zoltan2 DEFAULT_MSG Zoltan2_INCLUDE_DIR Zoltan2_LIBRARY)

MARK_AS_ADVANCED(Zoltan2_INCLUDE_DIR Zoltan2_LIBRARY)
