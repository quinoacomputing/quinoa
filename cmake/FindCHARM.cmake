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
  ${CYGWIN_INSTALL_PATH}/bin
)

if(CHARM_COMPILER)
  _GET_CHARMINC(HINTS_CHARMINC ${CHARM_COMPILER})
endif()

FIND_PATH(CHARM_INCLUDE_DIR NAMES charm.h
          HINTS ${HINTS_CHARMINC} ${CHARM_ROOT}/include
          PATH_SUFFIXES charm)

if(NOT BUILD_SHARED_LIBS)
  find_library(CHARM_CONV_CPLUS_LIBRARY NAMES libconv-cplus-y.a
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
  find_library(CHARM_CONV_UTIL_LIBRARY NAMES libconv-util.a
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
  find_library(CHARM_CONV_CORE_LIBRARY NAMES libconv-core.a
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
else()
  find_library(CHARM_CONV_CPLUS_LIBRARY NAMES conv-cplus-y
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
  find_library(CHARM_CONV_UTIL_LIBRARY NAMES conv-util
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
  find_library(CHARM_CONV_CORE_LIBRARY NAMES conv-core
               HINTS ${CHARM_ROOT}/lib $ENV{CHARM_ROOT}/lib)
endif()

set(CHARM_LIBRARIES ${CHARM_CONV_CPLUS_LIBRARY}
                    ${CHARM_CONV_UTIL_LIBRARY}
                    ${CHARM_CONV_CORE_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CHARM DEFAULT_MSG CHARM_COMPILER CHARM_INCLUDE_DIR CHARM_LIBRARIES CHARM_CONV_CPLUS_LIBRARY CHARM_CONV_UTIL_LIBRARY CHARM_CONV_CORE_LIBRARY)

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIR CHARM_LIBRARIES CHARM_CONV_CPLUS_LIBRARY CHARM_CONV_UTIL_LIBRARY CHARM_CONV_CORE_LIBRARY)

