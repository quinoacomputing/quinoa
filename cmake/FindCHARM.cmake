# Copyright (C) 2017 Los Alamos Natl. Lab
#
# This file was derived from FindGnuplot.cmake shipped with CMake 2.6.3.
#
# - this module looks for charmc
#
# Once done this will define
#
#  CHARM_FOUND - system has charmc
#  CHARM_COMPILER - the charmc executable
#
#  Usage:
#
#  set(CHARM_ROOT "/path/to/custom/charm") # prefer over system
#  find_package(HCHARM)
#  if(CHARM_FOUND)
#    #use CHARM_COMPILER
#  endif()

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
function(_GET_CHARMINC _OUT_INC _charmc)
  file(STRINGS ${_charmc} _contents REGEX "^CHARMINC=")
  if(_contents)
    string(REGEX REPLACE "^CHARMINC=\"(.*)\"" "\\1" ${_OUT_INC} "${_contents}")
    set(${_OUT_INC} ${${_OUT_INC}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "file ${_charmc} does not exist")
  endif()
endfunction()

# If already in cache, be silent
if (CHARM_INCLUDE_DIR AND CHARM_COMPILER)
  set (CHARM_FIND_QUIETLY TRUE)
endif()

INCLUDE(FindCygwin)

FIND_PROGRAM(CHARM_COMPILER
  NAMES 
  charmc
  PATHS
  ${CHARM_ROOT}/bin
  $ENV{CHARM_ROOT}/bin
  ${CYGWIN_INSTALL_PATH}/bin
)

if(CHARM_COMPILER)
  _GET_CHARMINC(HINTS_CHARMINC ${CHARM_COMPILER})
endif()

FIND_PATH(CHARM_INCLUDE_DIR NAMES charm.h
                            HINTS ${HINTS_CHARMINC}
                                  ${CHARM_ROOT}/include
                                  $ENV{CHARM_ROOT}/include
                            PATH_SUFFIXES charm)

# handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CHARM DEFAULT_MSG CHARM_COMPILER CHARM_INCLUDE_DIR)

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIR)

