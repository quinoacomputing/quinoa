################################################################################
#
# \file      cmake/FindNumDiff.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find NumDiff
# \date      Mon 27 Mar 2017 08:26:37 AM MDT
#
################################################################################

# Find NumDiff
#
#  NUMDIFF_FOUND - System has numdiff
#  NUMDIFF_EXECUTABLE - The numdiff executable
#
#  Usage:
#
#  set(NUMDIFF_ROOT "/path/to/custom/numdiff") # prefer over system
#  find_package(NumDiff)

if(NUMDIFF_EXECUTABLE)
  # Already in cache, be silent
  set (NUMDIFF_FIND_QUIETLY TRUE)
endif()

FIND_PROGRAM(NUMDIFF_EXECUTABLE NAMES numdiff
                                PATHS ${NUMDIFF_ROOT}/bin
                                      $ENV{NUMDIFF_ROOT}/bin
                                PATH_SUFFIXES bin)

# Handle the QUIETLY and REQUIRED arguments and set NUMDIFF_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NumDiff DEFAULT_MSG NUMDIFF_EXECUTABLE)

get_filename_component(${NUMDIFF_EXECUTABLE} NUMDIFF_EXECUTABLE ABSOLUTE)
MARK_AS_ADVANCED(NUMDIFF_EXECUTABLE)
