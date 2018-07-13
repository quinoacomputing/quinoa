################################################################################
#
# \file      cmake/FindLAPACKE.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

if(NOT BUILD_SHARED_LIBS)
  find_library(LAPACKE_LIBRARY NAMES liblapacke.a HINTS ${LAPACKE_ROOT}/lib
                                                        $ENV{LAPACKE_ROOT}/lib
               PATH_SUFFIXES lapack lapacke)
  # Static linking requires explicitly all dependent libraries
  find_library(LAPACK_LIBRARY NAMES liblapack.a HINTS ${LAPACKE_ROOT}/lib
                                                      $ENV{LAPACKE_ROOT}/lib
               PATH_SUFFIXES lapack lapacke)
  find_library(BLAS_LIBRARY NAMES libblas.a HINTS ${LAPACKE_ROOT}/lib
                                                  $ENV{LAPACKE_ROOT}/lib
               PATH_SUFFIXES lapack lapacke)
  # These two are also searched on multiarch lib paths. This requires knowing
  # the architecture and the major version of the compiler. The compiler version
  # is detected in src/CMakeLists.txt and the architecture is detected by
  # including GNUInstallDirs in cmake/TPLs.cmake.
  find_library(GFORTRAN_LIBRARY NAMES libgfortran.a
               HINTS ${LAPACKE_ROOT}/lib
                     /usr/lib/gcc/${CMAKE_LIBRARY_ARCHITECTURE}/${CMAKE_CXX_COMPILER_MAJOR}
                     $ENV{LAPACKE_ROOT}
                     $ENV{LAPACKE_ROOT}/lib
               PATH_SUFFIXES lapack lapacke)
  if(NOT ARCH MATCHES "ppc64")
    find_library(QUADMATH_LIBRARY NAMES libquadmath.a
                 HINTS ${LAPACKE_ROOT}/lib
                       /usr/lib/gcc/${CMAKE_LIBRARY_ARCHITECTURE}/${CMAKE_CXX_COMPILER_MAJOR}
                       $ENV{LAPACKE_ROOT}/lib
                 PATH_SUFFIXES lapack lapacke)
  endif()
else()
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
endif()

set(LAPACKE_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${GFORTRAN_LIBRARY} ${QUADMATH_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
if(NOT BUILD_SHARED_LIBS)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACKE_LIBRARY
    LAPACK_LIBRARY BLAS_LIBRARY GFORTRAN_LIBRARY QUADMATH_LIBRARY
    LAPACKE_INCLUDE_DIRS)
else()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACKE_LIBRARY
    LAPACK_LIBRARY BLAS_LIBRARY LAPACKE_INCLUDE_DIRS)
endif()

MARK_AS_ADVANCED(LAPACKE_INCLUDE_DIRS LAPACKE_LIBRARIES)
