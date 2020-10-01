################################################################################
#
# \file      FindExaM2M.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the ExaM2M mesh-to-mesh trasnfer library
#
################################################################################

#  ExaM2M: https://github.com/ThePhD/sol2
#
#  EXAM2M_FOUND - System has ExaM2M
#  EXAM2M_INCLUDE_DIRS - The ExaM2M include directory
#
#  Set the EXAM2M_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path
#
#  Usage:
#
#  set(EXAM2M_ROOT "/path/to/custom/exam2m") # prefer over system
#  find_package(ExaM2M)
#  include_directories(${EXAM2M_INCLUDE_DIRS})

# If already in cache, be silent
if(EXAM2M_INCLUDE_DIRS)
  set (EXAM2M_FIND_QUIETLY TRUE)
endif()

# stub to force "found" until we have something to find
# comment this to trigger building ExaM2M
# set(EXAM2M_INCLUDE_DIR "/tmp/exam2m/exam2m.hpp")

# find_path(EXAM2M_INCLUDE_DIR NAMES exam2m.hpp
#                              HINTS ${EXAM2M_ROOT}
#                                    $ENV{EXAM2M_ROOT}
#                              PATH_SUFFIXES include)

set(EXAM2M_INCLUDE_DIRS ${EXAM2M_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set EXAM2M_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ExaM2M REQUIRED_VARS EXAM2M_INCLUDE_DIRS)

MARK_AS_ADVANCED(EXAM2M_INCLUDE_DIRS)
