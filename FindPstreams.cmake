################################################################################
#
# \file      FindLPstreams.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Pstreams library
#
################################################################################

# Pstreams: http://pstreams.sourceforge.net/
#
#  PSTREAMS_FOUND - System has pstreams
#  PSTREAMS_INCLUDE_DIRS - The pstreams include directories
#
#  Set PSTREAMS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(PSTREAMS_ROOT "/path/to/custom/pstreams") # prefer over system
#  find_package(Pstreams)

# If already in cache, be silent
if(PSTREAMS_INCLUDE_DIRS)
  set (PSTREAMS_FIND_QUIETLY TRUE)
endif()

FIND_PATH(PSTREAMS_INCLUDE_DIR NAMES pstream.h
          PATH_SUFFIXES pstreams
          HINTS ${PSTREAMS_ROOT}/include $ENV{PSTREAMS_ROOT})

if(PSTREAMS_INCLUDE_DIR)
  set(PSTREAMS_INCLUDE_DIRS ${PSTREAMS_INCLUDE_DIR})
else()
  set(PSTREAMS_INCLUDE_DIRS "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set PSTREAMS_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pstreams DEFAULT_MSG PSTREAMS_INCLUDE_DIRS)

MARK_AS_ADVANCED(PSTREAMS_INCLUDE_DIRS)
