################################################################################
#
# \file      FindCBLAS.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the C-interface to BLAS
#
################################################################################

# Find the C-interface to BLAS
#
#  CBLAS_FOUND - System has CBLAS
#  CBLAS_INCLUDE_DIRS - The CBLAS include directories
#  CBLAS_LIBRARIES - The libraries needed to use CBLAS with dynamic linking
#
#  Set CBLAS_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage (assuming dynamic linking):
#
#  set(CBLAS_ROOT "/path/to/custom/cblas") # prefer over system
#  find_package(CBLAS)
#  if(CBLAS_FOUND)
#    target_link_libraries (TARGET ${CBLAS_LIBRARIES})
#  endif()
#
# In case BLAS is not installed in the default directory, set the CBLAS_ROOT
# variable to point to the root of BLAS, such that 'cblas.h' can be found in
# $CBLAS_ROOT/include. This can either be done using an environmental variable
# (e.g. export CBLAS_ROOT=/path/to/BLAS) or using a CMake variable
# (e.g. cmake -DCBLAS_ROOT=/path/to/BLAS ..).

# If already in cache, be silent
if(CBLAS_LIBRARIES AND CBLAS_INCLUDE_DIRS)
  set (CBLAS_FIND_QUIETLY TRUE)
endif()

find_path(CBLAS_INCLUDE_DIR cblas.h DOC "C-interface to BLAS"
          HINTS ${CBLAS_ROOT}/include $ENV{CBLAS_ROOT}/include
          PATH_SUFFIXES blas cblas)

find_library(CBLAS_LIBRARY NAMES cblas blas refblas
             HINTS ${CBLAS_ROOT}/lib
                   $ENV{CBLAS_ROOT}/lib
             PATH_SUFFIXES blas cblas)

set(CBLAS_INCLUDE_DIRS ${CBLAS_INCLUDE_DIR})
set(CBLAS_LIBRARIES ${CBLAS_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set CBLAS_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(CBLAS DEFAULT_MSG CBLAS_LIBRARY
  CBLAS_INCLUDE_DIRS)

MARK_AS_ADVANCED(CBLAS_INCLUDE_DIRS CBLAS_LIBRARIES)
