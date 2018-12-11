################################################################################
#
# \file      cmake/FindHypre.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the Hypre library from LLNL
#
################################################################################

# Find the Hypre library from LLNL
#
#  HYPRE_FOUND - System has Hypre
#  HYPRE_INCLUDE_DIRS - The Hypre include directory
#  HYPRE_LIBRARIES - The libraries needed to use Hypre
#
#  Set HYPRE_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(HYPRE_ROOT "/path/to/custom/hypre") # prefer over system
#  find_package(Hypre)
#  if(HYPRE_FOUND)
#    target_link_libraries (TARGET ${HYPRE_LIBRARIES})
#  endif()

function(_HYPRE_GET_VERSION _OUT_ver _version_hdr)
  file(STRINGS ${_version_hdr} _contents REGEX "#define HYPRE_RELEASE_VERSION[ \t]+")
  if(_contents)
    string(REGEX REPLACE "\"" "" _cont "${_contents}")
    string(REGEX REPLACE ".*#define HYPRE_RELEASE_VERSION[ \t]+([0-9.]+).*" "\\1" ${_OUT_ver} "${_cont}")
    if(NOT ${${_OUT_ver}} MATCHES "[0-9]+")
        message(FATAL_ERROR "Version parsing failed for HYPRE_RELEASE_VERSION in ${_version_hdr}!")
    endif()
    set(${_OUT_ver} ${${_OUT_ver}} PARENT_SCOPE)
 elseif(EXISTS ${_version_hdr})
    message(FATAL_ERROR "No HYPRE_RELEASE_VERSION in ${_version_hdr}")
 else()
    message(FATAL_ERROR "Include file ${_version_hdr} does not exist")
  endif()
endfunction()

# If already in cache, be silent
if(HYPRE_INCLUDE_DIRS AND HYPRE_LIBRARIES)
  set (HYPRE_FIND_QUIETLY TRUE)
endif()

if (HYPRE_ROOT)
  set(HYPRE_SEARCH_OPTS NO_DEFAULT_PATH)
endif()

find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h
                            PATH_SUFFIXES hypre
                            HINTS ${HYPRE_ROOT}/include
                            ${HYPRE_SEARCH_OPTS})

if(HYPRE_INCLUDE_DIR)
  _HYPRE_GET_VERSION(HYPRE_VERSION ${HYPRE_INCLUDE_DIR}/HYPRE_config.h)
  set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
else()
  set(HYPRE_VERSION 0.0.0)
  set(HYPRE_INCLUDE_DIRS "")
endif()

if(NOT BUILD_SHARED_LIBS)
  find_library(HYPRE_LIBRARY NAMES libHYPRE.a HINTS ${HYPRE_ROOT}/lib)
else()
  find_library(HYPRE_LIBRARY NAMES HYPRE HINTS ${HYPRE_ROOT}/lib)
endif()

set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set HYPRE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Hypre REQUIRED_VARS HYPRE_LIBRARIES HYPRE_INCLUDE_DIRS VERSION_VAR HYPRE_VERSION)

MARK_AS_ADVANCED(HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)
