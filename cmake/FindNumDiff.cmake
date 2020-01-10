################################################################################
#
# \file      FindNumDiff.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find NumDiff
#
################################################################################

# NumDiff: https://github.com/MethodicalAcceleratorDesign/MAD-X/tree/master/tools/numdiff
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

FIND_PROGRAM(NUMDIFF_EXECUTABLE NAMES maddiff mad-numdiff numdiff
                                PATHS ${NUMDIFF_ROOT}/bin
                                      $ENV{NUMDIFF_ROOT}/bin
                                PATH_SUFFIXES bin)

# Handle the QUIETLY and REQUIRED arguments and set NUMDIFF_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NumDiff DEFAULT_MSG NUMDIFF_EXECUTABLE)

get_filename_component(${NUMDIFF_EXECUTABLE} NUMDIFF_EXECUTABLE ABSOLUTE)
MARK_AS_ADVANCED(NUMDIFF_EXECUTABLE)
