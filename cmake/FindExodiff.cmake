################################################################################
#
# \file      FindExodiff.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find Exodiff
#
################################################################################

# Find the exodiff executable part of the ExodusII / Trilinos/SEACAS package
#
#  EXODIFF_FOUND      - True if the exodiff executable was found
#  EXODIFF_EXECUTABLE - The exodiff executable
#
#  Set EXODUS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(EXODUS_ROOT "/path/to/custom/exodus") # prefer over system
#  find_package(Exodiff)
#  if(EXODIFF_FOUND)
#    # use ${EXODIFF_EXECUTABLE} ...
#  endif()

INCLUDE(FindCygwin)

FIND_PROGRAM(EXODIFF_EXECUTABLE
  NAMES exodiff
  PATHS ${CYGWIN_INSTALL_PATH}/bin ${EXODUS_ROOT}/bin $ENV{EXODUS_ROOT}/bin
)

# handle the QUIETLY and REQUIRED arguments and set EXODIFF_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Exodiff DEFAULT_MSG EXODIFF_EXECUTABLE)

MARK_AS_ADVANCED(EXODIFF_EXECUTABLE)
