################################################################################
#
# \file      cmake/FindAEC.cmake
# \copyright 2016-2017, Los Alamos National Security, LLC.
# \brief     Find the Adaptive Entropy Coding library
#
################################################################################

# Find the Adaptive Entropy Coding library from
# https://www.dkrz.de/redmine/projects/aec/wiki
#
#  AEC_FOUND - System has AEC
#  AEC_INCLUDE_DIRS - The AEC include directory
#  AEC_LIBRARIES - The libraries needed to use AEC
#
#  Set AEC_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(AEC_ROOT "/path/to/custom/aec") # prefer over system
#  find_package(AEC)
#  if(AEC_FOUND)
#    target_link_libraries (TARGET ${AEC_LIBRARIES})
#  endif()

# If already in cache, be silent
if(AEC_INCLUDE_DIRS AND AEC_LIBRARY)
  set (AEC_FIND_QUIETLY TRUE)
endif()

find_path(AEC_INCLUDE_DIR NAMES szlib.h HINTS ${AEC_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
  find_library(AEC_LIBRARY NAMES libaec.a HINTS ${AEC_ROOT}/lib)
else()
  find_library(AEC_LIBRARY NAMES aec HINTS ${AEC_ROOT}/lib)
endif()

set(AEC_INCLUDE_DIRS ${AEC_INCLUDE_DIR})
set(AEC_LIBRARIES ${AEC_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set AEC_FOUND to TRUE if all
# listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(AEC DEFAULT_MSG AEC_LIBRARIES AEC_INCLUDE_DIRS)

MARK_AS_ADVANCED(AEC_INCLUDE_DIRS AEC_LIBRARIES)
