################################################################################
#
# \file      FindBackwardCpp.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the Backward-cpp header-only library
# \note      BackwardCPP requires other libraries, and for a static build,
# those libraries also require other libraries, e.g., elf, ebl, bz2, and lzma.
#
################################################################################

# Find Backward-cpp library headers
#
#  BACKWARDCPP_FOUND               - True if Backward-cpp is found
#  BACKWARD_INCLUDE_DIRS           - Backward-cpp include paths
#  BACKWARD_LIBRARIES              - Backward-cpp libraries to link
#
#  Set BACKWARDCPP_ROOT before calling find_package to a path to add an
#  additional search path, e.g.,
#
#  Usage:
#
#  set(BACKWARDCPP_ROOT "/path/to/custom/backward-cpp") # prefer over system
#  find_package(Backward)
#  if(BackwardCpp_FOUND)
#    include_directories( ${BACKWARD_INCLUDE_DIRS} )
#    target_link_libraries( <executable> ${BACKWARD_LIBRARIES} )
#  endif()

# If already in cache, be silent
if (BACKWARD_INCLUDE_DIRS AND BACKWARD_LIBRARIES)
  set (BACKWARDCPP_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(BACKWARDCPP_INCLUDE_DIR
          NAMES backward.hpp
          PATH_SUFFIXES include
          HINTS ${BACKWARDCPP_ROOT} $ENV{BACKWARDCPP_ROOT})

# Look for the cmake configuration file
FIND_PATH(BACKWARD_CMAKE_CONFIG_DIR
          NAMES BackwardConfig.cmake
          PATH_SUFFIXES lib/backward
          HINTS ${BACKWARDCPP_ROOT} $ENV{BACKWARDCPP_ROOT})

if(BACKWARDCPP_INCLUDE_DIR AND BACKWARD_CMAKE_CONFIG_DIR)
  list(APPEND CMAKE_MODULE_PATH "${BACKWARD_CMAKE_CONFIG_DIR}")
  include(BackwardConfig)

  # If BackwardCpp uses libdw, it needs additional libs for static builds
  if(NOT BUILD_SHARED_LIBS AND BACKWARD_LIBRARIES MATCHES "libdw")
    set(BACKWARD_STATIC_LIBS elf ebl bz2 lzma)
    foreach(lib ${BACKWARD_STATIC_LIBS})
      find_library(BACKWARD_${lib}_LIBRARY NAMES lib${lib}.a
                   HINTS /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE})
      list(APPEND BACKWARD_LIBRARIES "${BACKWARD_${lib}_LIBRARY}")
    endforeach()
  endif()

  message(STATUS "Backward-cpp config: ${BACKWARD_DEFINITIONS}")
  if (BACKWARD_LIBRARIES)
    message(STATUS "Backward-cpp libraries: ${BACKWARD_LIBRARIES}")
  endif()
endif()

set(BACKWARD_INCLUDE_DIRS ${BACKWARDCPP_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set BACKWARDCPP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BackwardCpp DEFAULT_MSG BACKWARD_INCLUDE_DIRS BACKWARD_LIBRARIES)
MARK_AS_ADVANCED(BACKWARD_INCLUDE_DIRS BACKWARD_LIBRARIES)
