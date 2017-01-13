################################################################################
#
# \file      cmake/FindRandom123.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find Random123
# \date      Wed 11 Jan 2017 12:15:51 PM MST
#
################################################################################

# Random123:
# http://www.thesalmons.org/john/random123/releases/latest/docs/index.html
#
#  Random123_FOUND - System has Random123
#  Random123_INCLUDE_PATH - The Random123 include directory
#
#  Set the RANDOM123_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(RANDOM123_ROOT "/path/to/custom/random123") # prefer over system
#  find_package(Random123)
#  include_directories(${RANDOM123_INCLUDE_PATH})

# If already in cache, be silent
if(RANDOM123_INCLUDE_PATH)
  set (RANDOM123_FIND_QUIETLY TRUE)
endif()

find_path(RANDOM123_INCLUDE_PATH NAMES philox.h
                            HINTS ${RANDOM123_ROOT}/include
                                  $ENV{RANDOM123_ROOT}
                            PATH_SUFFIXES Random123)

# Handle the QUIETLY and REQUIRED arguments and set Random123_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Random123 DEFAULT_MSG RANDOM123_INCLUDE_PATH)

MARK_AS_ADVANCED(RANDOM123_INCLUDE_PATH)
