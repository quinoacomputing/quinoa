################################################################################
#
# \file      FindMutationPP.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Mutation++ chemical reactions library
#
################################################################################
#
# Exports:
#   MutationPP_FOUND
#   MUTATIONPP_INCLUDE_DIRS
#   MUTATIONPP_LIBRARIES

# If already found, be quiet
if(MUTATIONPP_INCLUDE_DIRS AND MUTATIONPP_LIBRARIES AND TARGET MutationPP::MutationPP)
  set(MutationPP_FIND_QUIETLY TRUE)
endif()

# Collect hint roots (add more if your project uses different vars)
set(_MUTATIONPP_HINTS
  "${MUTATIONPP_ROOT}"
  "${MutationPP_ROOT}"
  "${MPP_ROOT}"
  "${TPL_DIR}"
  ${CMAKE_PREFIX_PATH}
)

# ---- Headers ---------------------------------------------------------------
# Mutation++ header layouts vary; try multiple plausible header names + suffixes.
find_path(MUTATIONPP_INCLUDE_DIR
  NAMES
    mutation++.h
    mutationpp.h
    Mutation++.h
  HINTS ${_MUTATIONPP_HINTS}
  PATH_SUFFIXES
    include
    include/mutation++
    include/mutationpp
    include/Mutationpp
    mutation++/include
    mutationpp/include
    Mutationpp/include
)

# Fallback: some installs put headers under .../include and require including <mutation++/mutation++.h>
# In that case MUTATIONPP_INCLUDE_DIR should be the parent include dir.
# If the header was found under an include subdir, CMake already returns the parent where the header was found.

# ---- Library ---------------------------------------------------------------
# Library names can differ across distros/builds.
find_library(MUTATIONPP_LIBRARY
  NAMES
    mutation++
    mutationpp
    mutation++_shared
    mutationpp_shared
    libmutation++
    libmutationpp
  HINTS ${_MUTATIONPP_HINTS}
  PATH_SUFFIXES
    lib
    lib64
    lib/static
    lib/shared
    mutation++/lib
    mutationpp/lib
    Mutationpp/lib
)

# ---- Optional data dir -----------------------------------------------------
unset(MUTATIONPP_DATA_DIR CACHE)
foreach(_h IN LISTS _MUTATIONPP_HINTS)
  if(NOT _h)
    continue()
  endif()
  foreach(_sfx
      share/mutation++/data
      share/mutationpp/data
      share/Mutationpp/data
      data
    )
    if(EXISTS "${_h}/${_sfx}")
      set(MUTATIONPP_DATA_DIR "${_h}/${_sfx}")
      break()
    endif()
  endforeach()
  if(MUTATIONPP_DATA_DIR)
    break()
  endif()
endforeach()

# ---- Standard args ---------------------------------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MutationPP
  REQUIRED_VARS MUTATIONPP_INCLUDE_DIR MUTATIONPP_LIBRARY
)

# ---- Export variables + imported target ------------------------------------
if(MutationPP_FOUND)
  set(MUTATIONPP_INCLUDE_DIRS "${MUTATIONPP_INCLUDE_DIR}")
  set(MUTATIONPP_LIBRARIES    "${MUTATIONPP_LIBRARY}")

  # Provide the target your project is already linking against.
  if(NOT TARGET MutationPP::MutationPP)
    add_library(MutationPP::MutationPP UNKNOWN IMPORTED)
    set_target_properties(MutationPP::MutationPP PROPERTIES
      IMPORTED_LOCATION "${MUTATIONPP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MUTATIONPP_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(MUTATIONPP_INCLUDE_DIR MUTATIONPP_LIBRARY MUTATIONPP_DATA_DIR)