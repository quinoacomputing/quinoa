# Find the C-interface to LAPACK as well as LAPACK/BLAS
#
#  LAPACKE_FOUND - System has LAPACKE as well as LAPACK/BLAS
#  LAPACKE_INCLUDE - The LAPACKE include directory
#  LAPACKE_LIBRARY - The libraries needed to use LAPACKE with dynamic linking
#  LAPACK_LIBRARIES - The libraries required to link LAPACKE with static linking
#
#  Set LAPACKE_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage (assuming dynamic linking):
#
#  set(LAPACKE_ROOT "/path/to/custom/lapacke") # prefer over system
#  find_package(LAPACKE)
#  if(LAPACKE_FOUND)
#    target_link_libraries (TARGET ${LAPACKE_LIBRARY})
#  endif()

# If already in cache, be silent
if(LAPACKE_INCLUDE AND LAPACKE_LIBRARY AND LAPACK_LIBRARIES)
  set (LAPACKE_FIND_QUIETLY TRUE)
endif()

find_path(LAPACKE_INCLUDE lapacke.h DOC "C-interface to LAPACK"
          HINTS ${LAPACKE_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
  find_library(LAPACKE_LIBRARY NAMES liblapacke.a HINTS ${LAPACKE_ROOT}/lib)
  # Static linking requires explicitly all dependent libraries
  find_library(LAPACK_LIBRARY NAMES liblapack.a HINTS ${LAPACKE_ROOT}/lib)
  find_library(BLAS_LIBRARY NAMES libblas.a HINTS ${LAPACKE_ROOT}/lib)
  # These two are also searched on multiarch lib paths. This requires knowing
  # the architecture and the major version of the compiler. The compiler version
  # is detected in src/CMakeLists.txt and the architecture is detected by
  # including GNUInstallDirs in cmake/TPLs.cmake.
  find_library(GFORTRAN_LIBRARY NAMES libgfortran.a HINTS ${LAPACKE_ROOT}/lib
               /usr/lib/gcc/${CMAKE_LIBRARY_ARCHITECTURE}/${GCC_MAJOR})
  find_library(QUADMATH_LIBRARY NAMES libquadmath.a  HINTS ${LAPACKE_ROOT}/lib
               /usr/lib/gcc/${CMAKE_LIBRARY_ARCHITECTURE}/${GCC_MAJOR})
else()
  find_library(LAPACKE_LIBRARY NAMES lapacke HINTS ${LAPACKE_ROOT}/lib)
endif()

# Only for reporting on what is found
set(LAPACK_LIBRARIES ${LAPACKE_LIBRARY} ${LAPACK_LIBRARY} ${BLAS_LIBRARY}
                     ${GFORTRAN_LIBRARY} ${QUADMATH_LIBRARY})

# For output only
set(LAPACK_LIBS ${LAPACK_LIBRARIES})

# Handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACK_LIBS)

MARK_AS_ADVANCED(LAPACKE_INCLUDE LAPACKE_LIBRARY LAPACK_LIBRARIES)
