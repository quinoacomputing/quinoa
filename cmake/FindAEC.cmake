################################################################################
#
# \file      cmake/FindAEC.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Adaptive Entropy Coding library
# \date      Fri 06 May 2016 06:41:51 AM MDT
#
################################################################################

# Find the Adaptive Entropy Coding library from
# https://www.dkrz.de/redmine/projects/aec/wiki
#
#  AEC_FOUND - System has AEC
#  AEC_INCLUDE - The AEC include directory
#  AEC_LIBRARY - The libraries needed to use AEC
#
#  Set AEC_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(AEC_ROOT "/path/to/custom/aec") # prefer over system
#  find_package(AEC)
#  if(AEC_FOUND)
#    target_link_libraries (TARGET ${AEC_LIBRARY})
#  endif()

# If already in cache, be silent
if(AEC_INCLUDE AND AEC_LIBRARY)
  set (AEC_FIND_QUIETLY TRUE)
endif()

find_path(AEC_INCLUDE NAMES szlib.h HINTS ${AEC_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
  find_library(AEC_LIBRARY NAMES libaec.a HINTS ${AEC_ROOT}/lib)
else()
  find_library(AEC_LIBRARY NAMES aec HINTS ${AEC_ROOT}/lib)
endif()

set(AEC_LIBRARY ${AEC_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set AEC_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(AEC DEFAULT_MSG AEC_LIBRARY AEC_INCLUDE)

MARK_AS_ADVANCED(AEC_INCLUDE AEC_LIBRARY)
