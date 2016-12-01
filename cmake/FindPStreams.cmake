################################################################################
#
# \file      cmake/FindLPstreams.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Pstreams library
# \date      Thu 01 Dec 2016 03:35:02 PM MST
#
################################################################################

# Find the Pstreams library
#
#  PSTREAMS_FOUND - System has pstreams
#  PSTREAMS_INCLUDE_PATH - The pstreams include directory (pstreams is header-only)
#
#  Set PSTREAMS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(PSTREAMS_ROOT "/path/to/custom/pstreams") # prefer over system
#  find_package(Pstreams)

# If already in cache, be silent
if(PSTREAMS_INCLUDE_PATH)
  set (PSTREAMS_FIND_QUIETLY TRUE)
endif()

FIND_PATH(PSTREAMS_INCLUDE_PATH NAMES pstream.h
          PATH_SUFFIXES pstreams
          HINTS ${PSTREAMS_ROOT}/include $ENV{PSTREAMS_ROOT})

# Handle the QUIETLY and REQUIRED arguments and set PSTREAMS_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pstreams DEFAULT_MSG PSTREAMS_INCLUDE_PATH)

MARK_AS_ADVANCED(PSTREAMS_INCLUDE_PATH)
