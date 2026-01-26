################################################################################
#
# \file      cmake/FindPugixml.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find the pugixml library
#
################################################################################

# pugixml: http://pugixml.org
#
#  PUGIXML_FOUND - System has pugixml
#  PUGIXML_INCLUDE_DIRS - The pugixml include directory (pugixml is header-only)
#  PUGIXML_LIBRARIES - The libraries needed to use pugixml
#
#  Set PUGIXML_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(PUGIXML_ROOT "/path/to/custom/pugixml") # prefer over system
#  find_package(Pugixml)
#  ifPUGIXML_FOUND)
#    target_link_libraries (TARGET ${PUGIXML_LIBRARIES})
#  endif()

# If already in cache, be silent
if(PUGIXML_INCLUDE_DIRS)
  set (PUGIXML_FIND_QUIETLY TRUE)
endif()

FIND_PATH(PUGIXML_INCLUDE_DIR NAMES pugixml.hpp HINTS ${PUGIXML_ROOT}/include
                                                $ENV{PUGIXML_ROOT}/include)

find_library(PUGIXML_LIBRARY NAMES pugixml HINTS ${PUGIXML_ROOT}
                                                 $ENV{PUGIXML_ROOT}
                                           PATH_SUFFIXES lib lib64)

if(PUGIXML_INCLUDE_DIR)
  set(PUGIXML_INCLUDE_DIRS ${PUGIXML_INCLUDE_DIR})
else()
  set(PUGIXML_INCLUDE_DIRS "")
endif()

set(PUGIXML_LIBRARIES ${PUGIXML_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set PUGIXML_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pugixml DEFAULT_MSG PUGIXML_LIBRARIES PUGIXML_INCLUDE_DIRS)

MARK_AS_ADVANCED(PUGIXML_INCLUDE_DIRS PUGIXML_LIBRARIES)
