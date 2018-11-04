################################################################################
#
# \file      cmake/FindHighwayHash.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find HighwayHash
#
################################################################################

# See HighwayHash: https://github.com/google/highwayhash
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

find_path(HIGHWAYHASH_INCLUDE_DIR NAMES sip_hash.h
                                  HINTS ${HIGHWAYHASH_ROOT}
                                        $ENV{HIGHWAYHASH_ROOT}
                                  PATH_SUFFIXES include highwayhash)

set(HIGHWAYHASH_INCLUDE_DIRS ${HIGHWAYHASH_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set HighwayHash_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HighwayHash REQUIRED_VARS HIGHWAYHASH_INCLUDE_DIRS)

MARK_AS_ADVANCED(HIGHWAYHASH_INCLUDE_DIRS)
