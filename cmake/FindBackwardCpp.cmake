################################################################################
#
# \file      cmake/FindBackwardCpp.cmake
# \copyright 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the Backward-cpp header-only library
# \note      BackwardCPP requires other libraries, and for a static build,
# those libraries also require other libraries, e.g., elf, ebl, bz2, and lzma.
#
################################################################################

# Find Backward-cpp library headers
#
#  BACKWARDCPP_FOUND               - True if Backward-cpp is found
#  BACKWARDCPP_INCLUDE_DIRS        - Backward-cpp install path
#  BACKWARDCPP_CMAKE_CONFIG_DIRS   - Backward-cpp cmake configuration path
#
#  Set BACKWARDCPP_ROOT before calling find_package to a path to add an
#  additional search path, e.g.,
#
#  Usage:
#
#  set(BACKWARDCPP_ROOT "/path/to/custom/backward-cpp") # prefer over system
#  find_package(Backward)
#  if(BackwardCpp_FOUND)
#    list(APPEND CMAKE_MODULE_PATH "${BACKWARDCPP_CMAKE_CONFIG_DIRS}")
#    include(BackwardConfig)
#    include_directories( ${BACKWARDCPP_INCLUDE_DIRS} )
#    target_link_libraries( <executable> ${BACKWARD_LIBRARIES} )
#  endif()

# If already in cache, be silent
if (BACKWARDCPP_INCLUDE_DIRS AND BACKWARDCPP_CMAKE_CONFIG_DIRS)
  set (BACKWARDCPP_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(BACKWARDCPP_INCLUDE_DIR
          NAMES backward.hpp
          PATH_SUFFIXES include
          HINTS ${BACKWARDCPP_ROOT} $ENV{BACKWARDCPP_ROOT})

# Look for the cmake configuration file
FIND_PATH(BACKWARDCPP_CMAKE_CONFIG_DIR
          NAMES BackwardConfig.cmake
          PATH_SUFFIXES lib/backward
          HINTS ${BACKWARDCPP_ROOT} $ENV{BACKWARDCPP_ROOT})

set(BACKWARDCPP_INCLUDE_DIRS ${BACKWARDCPP_INCLUDE_DIR})
set(BACKWARDCPP_CMAKE_CONFIG_DIRS ${BACKWARDCPP_CMAKE_CONFIG_DIR})

if(NOT BUILD_SHARED_LIBS AND BACKWARD_LIBRARIES MATCHES "libdw")
  set(BACKWARD_STATIC_LIBS elf ebl bz2 lzma)
  foreach(lib ${BACKWARD_STATIC_LIBS})
    find_library(BACKWARD_${lib}_LIBRARY NAMES lib${lib}.a
                 HINTS /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE})
    list(APPEND BACKWARD_STATIC_LIBRARIES "${BACKWARD_${lib}_LIBRARY}")
    list(APPEND BACKWARD_STATIC_LIBRARY_VARS BACKWARD_${lib}_LIBRARY)
  endforeach()
endif()

# Handle the QUIETLY and REQUIRED arguments and set BACKWARDCPP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
if(NOT BUILD_SHARED_LIBS)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(BackwardCpp DEFAULT_MSG
    BACKWARDCPP_INCLUDE_DIRS BACKWARDCPP_CMAKE_CONFIG_DIRS
    ${BACKWARD_STATIC_LIBRARY_VARS})
else()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(BackwardCpp DEFAULT_MSG
    BACKWARDCPP_INCLUDE_DIRS BACKWARDCPP_CMAKE_CONFIG_DIRS)
endif()

if (BACKWARD_STATIC_LIBRARIES)
  message(STATUS "Bacward-cpp static libraries: ${BACKWARD_STATIC_LIBRARIES}")
endif()

MARK_AS_ADVANCED(BACKWARDCPP_INCLUDE_DIRS BACKWARDCPP_CMAKE_CONFIG_DIRS
                 BACKWARD_STATIC_LIBRARIES)
