################################################################################
#
# \file      FindRandom123.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find Random123
#
################################################################################

# Random123:
# http://www.thesalmons.org/john/random123/releases/latest/docs/index.html
#
#  RANDOM123_FOUND - System has Random123
#  RANDOM123_INCLUDE_DIRS - The Random123 include directories
#
#  Set the RANDOM123_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(RANDOM123_ROOT "/path/to/custom/random123") # prefer over system
#  find_package(Random123)
#  include_directories(${RANDOM123_INCLUDE_DIRS})

# If already in cache, be silent
if(RANDOM123_INCLUDE_DIRS)
  set (RANDOM123_FIND_QUIETLY TRUE)
endif()

find_path(RANDOM123_INCLUDE_PHILOX NAMES Random123/philox.h
                                   HINTS ${RANDOM123_ROOT}
                                         $ENV{RANDOM123_ROOT})

find_path(RANDOM123_INCLUDE_THREEFRY NAMES Random123/threefry.h
                                     HINTS ${RANDOM123_ROOT}
                                           $ENV{RANDOM123_ROOT})

find_path(RANDOM123_INCLUDE_UNIFORM NAMES Random123/uniform.hpp
                                    HINTS ${RANDOM123_ROOT}
                                          $ENV{RANDOM123_ROOT})

set(RANDOM123_INCLUDE_DIRS ${RANDOM123_INCLUDE_PHILOX}
                           ${RANDOM123_INCLUDE_THREEFRY}
                           ${RANDOM123_INCLUDE_UNIFORM})

list(REMOVE_DUPLICATES RANDOM123_INCLUDE_DIRS)

# Handle the QUIETLY and REQUIRED arguments and set RANDOM123_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Random123 DEFAULT_MSG RANDOM123_INCLUDE_DIRS)

MARK_AS_ADVANCED(RANDOM123_INCLUDE_DIRS)
