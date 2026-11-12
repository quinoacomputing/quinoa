################################################################################
#
# \file      TPLsQHex.cmake
# \copyright 2026-     Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find third-party libraries required by QuinoaHex
#
################################################################################

# Add TPL_DIR to modules directory for TPLs that provide cmake FIND_PACKAGE code
SET(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

# Include support for multiarch path names
include(GNUInstallDirs)

#### TPLs we attempt to find on the system #####################################

message(STATUS "------------------------------------------")

#### BLAS/LAPACK library with LAPACKE C-interface
find_package(LAPACKE)
find_package(CBLAS)

#### Boost
set(BOOST_INCLUDEDIR ${TPL_DIR}/include) # prefer ours
find_package(Boost 1.56.0)
if(Boost_FOUND)
  message(STATUS "Boost at ${Boost_INCLUDE_DIR} (include)")
  include_directories(${Boost_INCLUDE_DIR})
endif()

#### AMReX
set(AMREX_ROOT ${TPL_DIR})
find_path(AMREX_INCLUDE_DIR NAMES AMReX.H
          HINTS ${AMREX_ROOT}/include
          NO_DEFAULT_PATH)
find_library(AMREX_LIBRARY NAMES amrex amrex_3d
             HINTS ${AMREX_ROOT}/lib ${AMREX_ROOT}/lib64
             NO_DEFAULT_PATH)
if (AMREX_INCLUDE_DIR AND AMREX_LIBRARY)
  set(AMREX_FOUND true)
  set(AMREX_INCLUDE_DIRS ${AMREX_INCLUDE_DIR})
  set(AMREX_LIBRARIES ${AMREX_LIBRARY})
endif()

#### Configure Backward-cpp
set(BACKWARD_ROOT ${TPL_DIR}) # prefer ours
find_package(BackwardCpp)
if(BACKWARDCPP_FOUND)
  set(HAS_BACKWARD true)  # will become compiler define
  message(STATUS "BackwardCpp enabled")
else()
  set(BACKWARD_INCLUDE_DIRS "")
  set(BACKWARD_LIBRARIES "")
endif()

#### Configure Brigand
set(BRIGAND_ROOT ${TPL_DIR}) # prefer ours
find_package(Brigand)

#### Configure Sol2
set(SOL2_ROOT ${TPL_DIR}) # prefer ours

find_package(Lua)
find_package(Sol2)
if (LUA_FOUND AND Sol2_FOUND)
  set(HAS_LUA true)  # will become compiler define
  message(STATUS "Lua enabled")
else()
  set(LUA_INCLUDE_DIR "")
endif()

message(STATUS "------------------------------------------")

# Function to print a list of missing library names
# Arguments:
#   'target' a string to use in the error message printed for which libraries are not found
#   'reqlibs' list of cmake variables in the form of "AMREX_FOUND", etc.
# Details: For each variable in 'reqlibs' if evaluates to false, trim the
# ending "_FOUND", convert to lower case and print an error message with the
# list of missing variables names. Intended to use after multiple find_package
# calls, passing all cmake variables named '*_FOUND' for all required libraries
# for a target.
function(PrintMissing target reqlibs)
  foreach(lib ${reqlibs})
    if(NOT ${lib})
      string(REPLACE "_FOUND" "" lib ${lib})
      string(TOLOWER ${lib} lib)
      list(APPEND missing "${lib}")
    endif()
  endforeach()
  string(REPLACE ";" ", " missing "${missing}")
  message(STATUS "Target '${target}' will NOT be configured, missing: ${missing}")
endfunction(PrintMissing)

# Enable inciterhex based on required TPLs found
if (AMREX_FOUND AND BRIGAND_FOUND AND LAPACKE_FOUND AND Boost_FOUND)
  set(INCITERHEX_EXECUTABLE inciterhex)
  set(ENABLE_INCITERHEX true CACHE BOOL "Enable ${INCITERHEX_EXECUTABLE}")
  if (NOT ENABLE_INCITERHEX)
    message(STATUS "Target '${INCITERHEX_EXECUTABLE}' disabled")
  endif()
else()
  PrintMissing(inciterhex "AMREX_FOUND;BRIGAND_FOUND;LAPACKE_FOUND;Boost_FOUND")
endif()
