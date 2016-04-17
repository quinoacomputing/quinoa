# Find the C-interface to LAPACK
#
#  LAPACKE_FOUND - System has LAPACKE
#  LAPACKE_INCLUDE - The LAPACKE include directory
#  LAPACKE_LIBRARY - The libraries needed to use LAPACKE
#
#  Set LAPACKE_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(LAPACKE_ROOT "/path/to/custom/lapacke") # prefer over system
#  find_package(LAPACKE)
#  if(LAPACKE_FOUND)
#    target_link_libraries (TARGET ${LAPACKE_LIBRARY})
#  endif()

# If already in cache, be silent
if(LAPACKE_INCLUDE AND LAPACKE_LIBRARY)
  set (LAPACKE_FIND_QUIETLY TRUE)
endif()

find_path(LAPACKE_INCLUDE lapacke.h DOC "C-interface to LAPACK"
          HINTS ${LAPACKE_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
  find_library(LAPACKE_LIBRARY NAMES liblapacke.a HINTS ${LAPACKE_ROOT}/lib)
else()
  find_library(LAPACKE_LIBRARY NAMES lapacke HINTS ${LAPACKE_ROOT}/lib)
endif()

# Handle the QUIETLY and REQUIRED arguments and set LAPACKE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LAPACKE DEFAULT_MSG LAPACKE_LIBRARY LAPACKE_INCLUDE)

MARK_AS_ADVANCED(LAPACKE_INCLUDE LAPACKE_LIBRARY)
