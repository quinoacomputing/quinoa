################################################################################
#
# \file      FindHighwayHash.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find HighwayHash
#
################################################################################

# HighwayHash: https://github.com/google/highwayhash
#
#  HighwayHash_FOUND - System has HighwayHash
#  HIGHWAYHASH_INCLUDE_DIRS - The HighwayHash include directory
#
#  Set the HIGHWAYHASH_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(HHIGHWAYHASH_ROOT "/path/to/custom/highwayhash") # prefer over system
#  find_package(HighwayHash)
#  include_directories(${HIGHWAYHASH_INCLUDE_DIRS})

# If already in cache, be silent
if(HIGHWAYHASH_INCLUDE_DIRS)
  set (HIGHWAYHASH_FIND_QUIETLY TRUE)
endif()

find_path(HIGHWAYHASH_INCLUDE_DIR NAMES highwayhash/sip_hash.h
                                  HINTS ${HIGHWAYHASH_ROOT}
                                        $ENV{HIGHWAYHASH_ROOT})

set(HIGHWAYHASH_INCLUDE_DIRS ${HIGHWAYHASH_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set HighwayHash_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HighwayHash REQUIRED_VARS HIGHWAYHASH_INCLUDE_DIRS)

MARK_AS_ADVANCED(HIGHWAYHASH_INCLUDE_DIRS)
