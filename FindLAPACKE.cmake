################################################################################
#
# \file      FindLAPACKE.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the C-interface to LAPACK as well as LAPACK/BLAS
#
################################################################################

# Find the C-interface to LAPACK as well as LAPACK/BLAS
#
#  LAPACKE_FOUND - System has LAPACKE as well as LAPACK/BLAS
#  LAPACKE_INCLUDE_DIRS - The LAPACKE include directories
#  LAPACKE_LIBRARIES - The libraries needed to use LAPACKE with dynamic linking
#
#  Set LAPACKE_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage (assuming dynamic linking):
#
#  set(LAPACKE_ROOT "/path/to/custom/lapacke") # prefer over system
#  find_package(LAPACKE)
#  if(LAPACKE_FOUND)
#    target_link_libraries (TARGET ${LAPACKE_LIBRARIES})
#  endif()

# If already in cache, be silent
if(LAPACKE_LIBRARIES AND LAPACKE_INCLUDE_DIRS)
  set (LAPACKE_FIND_QUIETLY TRUE)
endif()

find_path(LAPACKE_INCLUDE_DIR lapacke.h DOC "C-interface to LAPACK"
          HINTS ${LAPACKE_ROOT}/include $ENV{LAPACKE_ROOT}/include
          PATH_SUFFIXES lapack lapacke)

find_library(LAPACKE_LIBRARY NAMES lapacke reflapacke
             HINTS ${LAPACKE_ROOT}/lib
                   $ENV{LAPACKE_ROOT}/lib
             PATH_SUFFIXES lapack lapacke)
find_library(LAPACK_LIBRARY NAMES lapack reflapack
             HINTS ${LAPACKE_ROOT}/lib $ENV{LAPACKE_ROOT}/lib
             PATH_SUFFIXES lapack lapacke)
find_library(BLAS_LIBRARY NAMES blas refblas
             HINTS ${LAPACKE_ROOT}/lib $ENV{LAPACKE_ROOT}/lib
             PATH_SUFFIXES lapack lapacke)
set(GFORTRAN_LIBRARY "")
set(QUADMATH_LIBRARY "")

set(LAPACKE_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${GFORTRAN_LIBRARY} ${QUADMATH_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACKE_LIBRARY
  LAPACK_LIBRARY BLAS_LIBRARY LAPACKE_INCLUDE_DIRS)

MARK_AS_ADVANCED(LAPACKE_INCLUDE_DIRS LAPACKE_LIBRARIES)
