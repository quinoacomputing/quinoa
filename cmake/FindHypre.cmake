################################################################################
#
# \file      cmake/FindHypre.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Hypre library from LLNL
# \date      Wed 30 Nov 2016 08:05:15 AM MST
#
################################################################################

# Find the Hypre library from LLNL
#
#  HYPRE_FOUND - System has Hypre
#  HYPRE_INCLUDE_PATH - The Hypre include directory
#  HYPRE_LIBRARY - The libraries needed to use Hypre
#
#  Set HYPRE_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(HYPRE_ROOT "/path/to/custom/hypre") # prefer over system
#  find_package(Hypre)
#  if(HYPRE_FOUND)
#    target_link_libraries (TARGET ${HYPRE_LIBRARY})
#  endif()

# If already in cache, be silent
if(HYPRE_INCLUDE_PATH AND HYPRE_LIBRARY)
  set (HYPRE_FIND_QUIETLY TRUE)
endif()

find_path(HYPRE_INCLUDE_PATH NAMES HYPRE.h
                             HINTS ${HYPRE_ROOT}/include
                             HINTS /usr/include/hypre)

if(NOT BUILD_SHARED_LIBS)
  find_library(HYPRE_LIBRARY NAMES libHYPRE.a HINTS ${HYPRE_ROOT}/lib)
else()
  find_library(HYPRE_LIBRARY NAMES HYPRE HINTS ${HYPRE_ROOT}/lib)
endif()

# Handle the QUIETLY and REQUIRED arguments and set HYPRE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Hypre DEFAULT_MSG HYPRE_LIBRARY
                                  HYPRE_INCLUDE_PATH)

MARK_AS_ADVANCED(HYPRE_INCLUDE_PATH HYPRE_LIBRARY)
