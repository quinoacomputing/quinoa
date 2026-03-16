################################################################################
#
# \file      TPLs.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find third-party libraries required
#
################################################################################

# Add TPL_DIR to modules directory for TPLs that provide cmake FIND_PACKAGE
# code, such as Trilinos
SET(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

# Include support for multiarch path names
include(GNUInstallDirs)

#### TPLs we attempt to find on the system #####################################

message(STATUS "------------------------------------------")

#### Charm++
if (NOT CHARM_ROOT)
  set(CHARM_ROOT ${TPL_DIR}/charm)
endif()
find_package(Charm)

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

#### TUT
find_package(TUT)

#### PugiXML
set(PUGIXML_ROOT ${TPL_DIR}) # prefer ours
find_package(Pugixml)

#### PEGTL
find_package(PEGTL 2.0.0)

### NetCDF
find_package(NetCDF)

#### Zoltan2 library
find_package(Zoltan2)

#### NumDiff executable
find_package(NumDiff)

#### ExodusII library
find_package(SEACASExodus)
set(EXODUS_ROOT ${TPL_DIR}) # prefer ours
find_package(Exodiff)

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

#### Configure HighwayHash
set(HIGHWAYHASH_ROOT ${TPL_DIR}) # prefer ours
find_package(HighwayHash)

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

# For Mutation++
include_directories(SYSTEM ${TPL_DIR}/include/eigen3)

find_package(MutationPP QUIET)
if (MutationPP_FOUND)
  message(STATUS "Found Mutation++: ${MUTATIONPP_LIBRARY}")
  # Treat Mutation++ headers as system headers so compiler warnings from them are supressed
  include_directories(SYSTEM ${MUTATIONPP_INCLUDE_DIRS})
else()
  message(STATUS "Mutation++ not found; chemistry features will be disabled")
endif()

message(STATUS "------------------------------------------")

# Function to print a list of missing library names
# Arguments:
#   'target' a string to use in the error message printed for which libraries are not found
#   'reqlibs' list of cmake variables in the form of "CHARM_FOUND", Boost_FOUND, etc.
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

# Enable individual executables based on required TPLs found

if (CHARM_FOUND AND PUGIXML_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND
    BRIGAND_FOUND AND TUT_FOUND AND PEGTL_FOUND AND Boost_FOUND AND
    HIGHWAYHASH_FOUND AND LAPACKE_FOUND
    AND ENABLE_TESTS)
  set(UNITTEST_EXECUTABLE unittest)
  set(ENABLE_UNITTEST true CACHE BOOL "Enable ${UNITTEST_EXECUTABLE}")
  if (NOT ENABLE_UNITTEST)
    message(STATUS "Target '${UNITTEST_EXECUTABLE}' disabled")
  endif()
else()
  if (NOT ENABLE_TESTS)
    message(STATUS "Target 'unittest' will NOT be configured, tests disabled.")
  else()
    PrintMissing(unittest "CHARM_FOUND;PUGIXML_FOUND;SEACASExodus_FOUND;EXODIFF_FOUND;BRIGAND_FOUND;TUT_FOUND;PEGTL_FOUND;Boost_FOUND;HIGHWAYHASH_FOUND;LAPACKE_FOUND")
  endif()
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND
    Zoltan2_FOUND AND BRIGAND_FOUND AND PEGTL_FOUND AND
    LAPACKE_FOUND AND Boost_FOUND AND HIGHWAYHASH_FOUND)
  set(INCITER_EXECUTABLE inciter)
  set(ENABLE_INCITER true CACHE BOOL "Enable ${INCITER_EXECUTABLE}")
  if (NOT ENABLE_INCITER)
    message(STATUS "Target '${INCITER_EXECUTABLE}' disabled")
  endif()
else()
  PrintMissing(inciter "CHARM_FOUND;SEACASExodus_FOUND;EXODIFF_FOUND;Zoltan2_FOUND;BRIGAND_FOUND;PEGTL_FOUND;LAPACKE_FOUND;Boost_FOUND")
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND PEGTL_FOUND AND
    PUGIXML_FOUND AND Boost_FOUND AND BRIGAND_FOUND AND
    HIGHWAYHASH_FOUND)
  set(MESHCONV_EXECUTABLE meshconv)
  set(ENABLE_MESHCONV true CACHE BOOL "Enable ${MESHCONV_EXECUTABLE}")
  if (NOT ENABLE_MESHCONV)
    message(STATUS "Target '${MESHCONV_EXECUTABLE}' disabled")
  endif()
else()
  PrintMissing(meshconv "CHARM_FOUND;SEACASExodus_FOUND;EXODIFF_FOUND;PEGTL_FOUND;PUGIXML_FOUND;Boost_FOUND;BRIGAND_FOUND;HIGHWAYHASH_FOUND")
endif()
