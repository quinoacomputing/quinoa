################################################################################
#
# \file      cmake/Findpugixml.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the pugixml library
# \date      Fri 13 Jan 2017 09:20:02 AM MST
#
################################################################################

# Find the pugixml library
#
#  PUGIXML_FOUND - System has pugixml
#  PUGIXML_INCLUDES - The pugixml include directory (pugixml is header-only)
#  PUGIXML_LIBRARIES - The libraries needed to use pugixml
#
#  Set PUGIXML_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(PUGIXML_ROOT "/path/to/custom/pugixml") # prefer over system
#  find_package(pugixml)
#  ifPUGIXML_FOUND)
#    target_link_libraries (TARGET ${PUGIXML_LIBRARIES})
#  endif()

# If already in cache, be silent
if(PUGIXML_INCLUDES)
  set (PUGIXML_FIND_QUIETLY TRUE)
endif()

FIND_PATH(PUGIXML_INCLUDES NAMES pugixml.hpp HINTS ${PUGIXML_ROOT}/include
                                                   $ENV{PUGIXML_ROOT}/include)

if(NOT BUILD_SHARED_LIBS)
 find_library(PUGIXML_LIBRARIES NAMES libpugixml.a HINTS ${PUGIXML_ROOT}/lib
                                                   $ENV{PUGIXML_ROOT}/lib64)
else()
 find_library(PUGIXML_LIBRARIES NAMES pugixml HINTS ${PUGIXML_ROOT}/lib
                                                    $ENV{PUGIXML_ROOT}/lib64)
endif()

# Handle the QUIETLY and REQUIRED arguments and set PUGIXML_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(pugixml DEFAULT_MSG PUGIXML_LIBRARIES PUGIXML_INCLUDES)

MARK_AS_ADVANCED(PUGIXML_INCLUDES PUGIXML_LIBRARIES)
