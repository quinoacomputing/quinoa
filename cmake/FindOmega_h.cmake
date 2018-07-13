################################################################################
#
# \file      cmake/FindOmega_h.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the Omega_h library
#
################################################################################

# Find Omega_h headers and libraries
#
#  OMEGA_H_FOUND        - True if Omega_h is found
#  OMEGA_H_INCLUDE_DIRS - Omega_h include files paths
#  OMEGA_H_LIBRARIES    - List of Omega_h libraries
#
#  Set OMEGA_H_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(OMEGA_H_ROOT "/path/to/custom/omega_h") # prefer over system
#  find_package(OMEGA_H)
#  if(OMEGA_H_FOUND)
#    target_link_libraries (TARGET ${OMEGA_H_LIBRARIES})
#  endif()

# If already in cache, be silent
if (OMEGA_H_INCLUDE_DIRS AND OMEGA_H_LIBRARIES)
  set (Omega_h_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(OMEGA_H_INCLUDE_DIR NAMES Omega_h_library.hpp
                              HINTS ${OMEGA_H_ROOT}/include
                                    $ENV{OMEGA_H_ROOT}/include)

# Look for the library
if(NOT BUILD_SHARED_LIBS)
  FIND_LIBRARY(OMEGA_H_LIBRARY NAMES libomega_h.a
                               HINTS ${OMEGA_H_ROOT}/lib
                                     $ENV{OMEGA_H_ROOT}/lib)
else()
  FIND_LIBRARY(OMEGA_H_LIBRARY NAMES omega_h
                               HINTS ${OMEGA_H_ROOT}/lib
                                     $ENV{OMEGA_H_ROOT}/lib)
endif()

set(OMEGA_H_INCLUDE_DIRS ${OMEGA_H_INCLUDE_DIR})
set(OMEGA_H_LIBRARIES ${OMEGA_H_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set OMEGA_H_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Omega_h DEFAULT_MSG OMEGA_H_INCLUDE_DIRS OMEGA_H_LIBRARIES)

MARK_AS_ADVANCED(OMEGA_H_INCLUDE_DIRS OMEGA_H_LIBRARIES)
