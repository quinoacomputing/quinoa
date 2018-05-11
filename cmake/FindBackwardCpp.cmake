################################################################################
#
# \file      cmake/FindBackwardCpp.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the Backward-cpp header-only library
#
################################################################################

# Find Backward-cpp library headers
#
#  BACKWARDCPP_FOUND               - True if Backward-cpp is found
#  BACKWARDCPP_INCLUDE_DIRS        - Backward-cpp install path
#  BACKWARDCPP_CMAKE_CONFIG_DIRS   - Backward-cpp cmake configuration path
#
#  Set BACKWARD_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(BACKWARD_ROOT "/path/to/custom/backward-cpp") # prefer over system
#  find_package(Backward)
#  if(BackwardCpp_FOUND)
#    list(APPEND CMAKE_MODULE_PATH "${BACKWARDCPP_CMAKE_CONFIG_DIRS}")
#    include(BackwardConfig)
#    include_directories( ${BACKWARDCPP_INCLUDE_DIRS} )
#    target_link_libraries( <executable> ${BACKWARD_LIBRARIES} )
#  endif()

# If already in cache, be silent
if (BACKWARDCPP_INCLUDE_DIRS AND BACKWARDCPP_CMAKE_CONFIG_DIRS)
  set (BACKWARD_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(BACKWARD_INCLUDE_DIR
          NAMES backward.hpp
          PATH_SUFFIXES include
          HINTS ${BACKWARD_ROOT} $ENV{BACKWARD_ROOT})

# Look for the cmake configuration file
FIND_PATH(BACKWARD_CMAKE_CONFIG_DIR
          NAMES BackwardConfig.cmake
          PATH_SUFFIXES lib/backward
          HINTS ${BACKWARD_ROOT} $ENV{BACKWARD_ROOT})

set(BACKWARDCPP_INCLUDE_DIRS ${BACKWARD_INCLUDE_DIR})
set(BACKWARDCPP_CMAKE_CONFIG_DIRS ${BACKWARD_CMAKE_CONFIG_DIR})

# Handle the QUIETLY and REQUIRED arguments and set BACKWARD_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Backward DEFAULT_MSG
  BACKWARDCPP_INCLUDE_DIRS BACKWARDCPP_CMAKE_CONFIG_DIRS)

MARK_AS_ADVANCED(BACKWARDCPP_INCLUDE_DIRS BACKWARDCPP_CMAKE_CONFIG_DIRS)
