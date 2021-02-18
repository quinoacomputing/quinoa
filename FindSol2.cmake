################################################################################
#
# \file      FindSol2.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Sol2 Lua C++ binding library
#
################################################################################

# Sol2: https://github.com/ThePhD/sol2
#
#  SOL2_FOUND - System has Sol2
#  SOL2_INCLUDE_DIRS - The Sol2 include directory
#
#  Set the SOL2_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(SOL2_ROOT "/path/to/custom/brigand") # prefer over system
#  find_package(Sol2)
#  include_directories(${SOL2_INCLUDE_DIRS})

# If already in cache, be silent
if(SOL2_INCLUDE_DIRS)
  set (SOL2_FIND_QUIETLY TRUE)
endif()

find_path(SOL2_INCLUDE_DIR NAMES sol.hpp
                           HINTS ${SOL2_ROOT}
                                 $ENV{SOL2_ROOT}
                           PATH_SUFFIXES sol include)

# Remove last 'sol' from path found, otherwise some compilers will not find
# sol.hpp.
get_filename_component(SOL2_INCLUDE_DIR ${SOL2_INCLUDE_DIR} DIRECTORY)

set(SOL2_INCLUDE_DIRS ${SOL2_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set SOL2_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Sol2 REQUIRED_VARS SOL2_INCLUDE_DIRS)

MARK_AS_ADVANCED(SOL2_INCLUDE_DIRS)
