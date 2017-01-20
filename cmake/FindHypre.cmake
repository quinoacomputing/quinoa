################################################################################
#
# \file      cmake/FindHypre.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Hypre library from LLNL
# \date      Fri 20 Jan 2017 12:19:59 PM MST
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

# If already in cache, be silent
if(HYPRE_INCLUDE_DIRS AND HYPRE_LIBRARIES)
  set (HYPRE_FIND_QUIETLY TRUE)
endif()

find_path(HYPRE_INCLUDE_DIR NAMES HYPRE.h
                            PATH_SUFFIXES hypre
                            HINTS ${HYPRE_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
  find_library(HYPRE_LIBRARY NAMES libHYPRE.a HINTS ${HYPRE_ROOT}/lib)
else()
  find_library(HYPRE_LIBRARY NAMES HYPRE HINTS ${HYPRE_ROOT}/lib)
endif()

set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set HYPRE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Hypre DEFAULT_MSG HYPRE_LIBRARIES
                                  HYPRE_INCLUDE_DIRS)

MARK_AS_ADVANCED(HYPRE_INCLUDE_DIRS HYPRE_LIBRARIES)
