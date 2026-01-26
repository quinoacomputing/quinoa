################################################################################
#
# \file      FindBrigand.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Brigand template metaprogramming library
#
################################################################################

# Brigand: https://github.com/edouarda/brigand
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

# Remove last 'brigand' from path found, otherwise some compilers, e.g., gc++
# 4.8, will not find brigand.
get_filename_component(BRIGAND_INCLUDE_DIR ${BRIGAND_INCLUDE_DIR} DIRECTORY)

set(BRIGAND_INCLUDE_DIRS ${BRIGAND_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set Brigand_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Brigand REQUIRED_VARS BRIGAND_INCLUDE_DIRS)

MARK_AS_ADVANCED(BRIGAND_INCLUDE_DIRS)
