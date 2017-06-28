################################################################################
#
# \file      cmake/FindH5Part.cmake
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the H5Part library
#
################################################################################

# Find H5Part headers and libraries
#
#  H5PART_FOUND        - True if H5Part is found
#  H5PART_INCLUDE_DIRS - H5Part include files paths
#  H5PART_LIBRARIES    - List of H5Part libraries
#
#  Set H5PART_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(H5PART_ROOT "/path/to/custom/h5part") # prefer over system
#  find_package(H5PART)
#  if(H5PART_FOUND)
#    target_link_libraries (TARGET ${H5PART_LIBRARIES})
#  endif()

# If already in cache, be silent
if (H5PART_INCLUDE_DIRS AND H5PART_LIBRARIES)
  set (H5Part_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(H5PART_INCLUDE_DIR NAMES H5Part.h HINTS ${H5PART_ROOT}/include
                                                  $ENV{H5PART_ROOT}/include)

# Look for the library
if(NOT BUILD_SHARED_LIBS)
  FIND_LIBRARY(H5PART_LIBRARY NAMES libH5Part.a HINTS ${H5PART_ROOT}/lib
                                                      $ENV{H5PART_ROOT}/lib)
else()
  FIND_LIBRARY(H5PART_LIBRARY NAMES H5Part HINTS ${H5PART_ROOT}/lib
                                                 $ENV{H5PART_ROOT}/lib)
endif()

set(H5PART_INCLUDE_DIRS ${H5PART_INCLUDE_DIR})
set(H5PART_LIBRARIES ${H5PART_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set H5PART_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(H5Part DEFAULT_MSG H5PART_INCLUDE_DIRS H5PART_LIBRARIES)

MARK_AS_ADVANCED(H5PART_INCLUDE_DIRS H5PART_LIBRARIES)
