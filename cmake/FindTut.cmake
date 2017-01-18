################################################################################
#
# \file      cmake/FindTut.cmake
# \author    C. Junghans
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Tut header
# \date      Wed 18 Jan 2017 04:11:59 PM MST
#
################################################################################

# Find Tut headers and libraries
#
#  Tut_FOUND        - True if Tut is found
#  Tut_INCLUDE_DIRS - Tut include files paths
#
#  Set Tut_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(Tut_ROOT "/path/to/custom/tut") # prefer over system
#  find_package(Tut)
#  if(Tut_FOUND)
#    include_directories(${Tut_INCLUDE_DIR})
#  endif()

find_package(PkgConfig)

# If already in cache, be silent
if (Tut_INCLUDE_DIR)
  set (Tut_FIND_QUIETLY TRUE)
  pkg_check_modules(PC_TUT QUIET tut)
else()
  pkg_check_modules(PC_TUT tut)
endif()

# Look for the header file
FIND_PATH(Tut_INCLUDE_DIR NAMES tut.h tut.hpp
                          HINTS ${Tut_ROOT}/include
                                $ENV{Tut_ROOT}/include 
                                ${PC_TUT_INCLUDE_DIRS}
                          PATH_SUFFIXES tut)

set(Tut_INCLUDE_DIRS ${Tut_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set Tut_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Tut DEFAULT_MSG Tut_INCLUDE_DIR)

MARK_AS_ADVANCED(Tut_INCLUDE_DIR)
