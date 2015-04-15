# Find the Hypre library from LLNL
#
#  HYPRE_FOUND - System has Hypre
#  HYPRE_INCLUDES - The Hypre include directory
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

if(HYPRE_INCLUDES AND HYPRE_LIBRARIES)
  # Already in cache, be silent
  set (HYPRE_FIND_QUIETLY TRUE)
endif()

FIND_PATH(HYPRE_INCLUDES NAMES HYPRE.h HINTS ${HYPRE_ROOT}/include)
FIND_LIBRARY(HYPRE_LIBRARIES NAMES HYPRE HINTS ${HYPRE_ROOT}/lib)

# Handle the QUIETLY and REQUIRED arguments and set HYPRE_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Hypre DEFAULT_MSG HYPRE_LIBRARIES HYPRE_INCLUDES)

MARK_AS_ADVANCED(HYPRE_INCLUDES HYPRE_LIBRARIES)
