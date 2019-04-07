################################################################################
#
# \file      FindPEGTL.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find PEGTL
#
################################################################################

# PEGTL: https://github.com/taocpp/PEGTL
#
#  PEGTL_FOUND - System has PEGTL
#  PEGTL_INCLUDE_DIRS - The PEGTL include directory
#
#  Set the PEGTL_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(PEGTL_ROOT "/path/to/custom/pegtl") # prefer over system
#  find_package(PEGTL)
#  include_directories(${PEGTL_INCLUDE_DIRS})

function(_PEGTL_GET_VERSION _OUT_major _OUT_minor _OUT_micro _metisversion_hdr)
    file(STRINGS ${_metisversion_hdr} _contents REGEX "#define TAOCPP_PEGTL_VERSION_[A-Z]+[ \t]+")
    if(_contents)
        string(REGEX REPLACE ".*#define TAOCPP_PEGTL_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" ${_OUT_major} "${_contents}")
        string(REGEX REPLACE ".*#define TAOCPP_PEGTL_VERSION_MINOR[ \t]+([0-9]+).*" "\\1" ${_OUT_minor} "${_contents}")
        string(REGEX REPLACE ".*#define TAOCPP_PEGTL_VERSION_PATCH[ \t]+([0-9]+).*" "\\1" ${_OUT_micro} "${_contents}")

        if(NOT ${_OUT_major} MATCHES "[0-9]+")
            message(FATAL_ERROR "Version parsing failed for TAOCPP_PEGTL_VERSION_MAJOR!")
        endif()
        if(NOT ${_OUT_minor} MATCHES "[0-9]+")
            message(FATAL_ERROR "Version parsing failed for TAOCPP_PEGTL_VERSION_MINOR!")
        endif()
        if(NOT ${_OUT_micro} MATCHES "[0-9]+")
            message(FATAL_ERROR "Version parsing failed for TAOCPP_PEGTL_VERSION_PATCH!")
        endif()

        set(${_OUT_major} ${${_OUT_major}} PARENT_SCOPE)
        set(${_OUT_minor} ${${_OUT_minor}} PARENT_SCOPE)
        set(${_OUT_micro} ${${_OUT_micro}} PARENT_SCOPE)

    else()
        message(FATAL_ERROR "Include file ${_metisversion_hdr} does not exist")
    endif()
endfunction()

# If already in cache, be silent
if(PEGTL_INCLUDE_DIRS)
  set (PEGTL_FIND_QUIETLY TRUE)
endif()

find_path(PEGTL_INCLUDE_DIR NAMES pegtl.hpp
                            HINTS ${PEGTL_ROOT}/include
                                  $ENV{PEGTL_ROOT}/include
                            PATH_SUFFIXES tao pegtl/include/tao)
if(PEGTL_INCLUDE_DIR)
  _PEGTL_GET_VERSION(PEGTL_MAJOR_VERSION PEGTL_MINOR_VERSION PEGTL_PATCH_VERSION ${PEGTL_INCLUDE_DIR}/pegtl/version.hpp)
  set(PEGTL_VERSION ${PEGTL_MAJOR_VERSION}.${PEGTL_MINOR_VERSION}.${PEGTL_PATCH_VERSION})
else()
  set(PEGTL_VERSION 0.0.0)
endif()

if(PEGTL_INCLUDE_DIR)
  set(PEGTL_INCLUDE_DIRS ${PEGTL_INCLUDE_DIR})
else()
  set(PEGTL_INCLUDE_DIRS "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set PEGTL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PEGTL REQUIRED_VARS PEGTL_INCLUDE_DIRS VERSION_VAR PEGTL_VERSION)

MARK_AS_ADVANCED(PEGTL_INCLUDE_DIRS)
