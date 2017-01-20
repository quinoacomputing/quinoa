################################################################################
#
# \file      cmake/FindExodiff.cmake
# \author    J. Bakosi
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Adaptive Entropy Coding library
# \date      Fri 20 Jan 2017 12:01:56 PM MST
#
################################################################################

# Find the exodiff executable part of the ExodusII / Trilinos/SEACAS package
#
#  Exodiff_FOUND      - True if the exodiff executable was found
#  EXODIFF_EXECUTABLE - The exodiff executable
#
#  Set EXODUS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(EXODUS_ROOT "/path/to/custom/exodus") # prefer over system
#  find_package(Exodiff)
#  if(Exodiff_FOUND)
#    # use ${EXODIFF_EXECUTABLE} ...
#  endif()

INCLUDE(FindCygwin)

FIND_PROGRAM(EXODIFF_EXECUTABLE
  NAMES exodiff
  PATHS ${CYGWIN_INSTALL_PATH}/bin ${EXODUS_ROOT}/bin $ENV{EXODUS_ROOT}/bin
)

# handle the QUIETLY and REQUIRED arguments and set Exodiff_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Exodiff DEFAULT_MSG EXODIFF_EXECUTABLE)

MARK_AS_ADVANCED(EXODIFF_EXECUTABLE)
