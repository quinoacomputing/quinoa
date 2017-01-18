# Copyright (C) 2016 Los Alamos Natl. Lab
#
# This file was derived from FindGnuplot.cmake shipped with CMake 2.6.3.
#
# - this module looks for exodiff
#
# Once done this will define
#
#  SEACASExodiff_FOUND - system has exodiff
#  SEACASExodiff_EXECUTABLE - the exodiff executable
#
#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

INCLUDE(FindCygwin)

FIND_PROGRAM(SEACASExodiff_EXECUTABLE
  NAMES exodiff
  PATHS ${CYGWIN_INSTALL_PATH}/bin $ENV{EXODUS_ROOT}/bin
)

# handle the QUIETLY and REQUIRED arguments and set SEACASExodiff_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SEACASExodiff DEFAULT_MSG SEACASExodiff_EXECUTABLE)

IF(NOT SEACASExodiff_FOUND)
  message(STATUS "exodiff not found, help cmake to find it by setting SEACASExodiff_EXECUTABLE")
ENDIF(NOT SEACASExodiff_FOUND)

MARK_AS_ADVANCED( SEACASExodiff_EXECUTABLE )

