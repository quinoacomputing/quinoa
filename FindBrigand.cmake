################################################################################
#
# \file      cmake/FindBrigand.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find Brigand
#
################################################################################

# See Brigand: https://github.com/edouarda/brigand
#
#  BRIGAND_FOUND - System has Brigand
#  BRIGAND_INCLUDE_DIRS - The Brigand include directory
#
#  Set the BRIGAND_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(BRIGAND_ROOT "/path/to/custom/brigand") # prefer over system
#  find_package(Brigand)
#  include_directories(${BRIGAND_INCLUDE_DIRS})

# If already in cache, be silent
if(BRIGAND_INCLUDE_DIRS)
  set (BRIGAND_FIND_QUIETLY TRUE)
endif()

find_path(BRIGAND_INCLUDE_DIR NAMES brigand.hpp
                              HINTS ${BRIGAND_ROOT}/include
                                    $ENV{BRIGAND_ROOT}/include
                              PATH_SUFFIXES brigand)

set(BRIGAND_INCLUDE_DIRS ${BRIGAND_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set Brigand_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Brigand REQUIRED_VARS BRIGAND_INCLUDE_DIRS)

MARK_AS_ADVANCED(BRIGAND_INCLUDE_DIRS)
