# Find the pugixml library
#
#  PUGIXML_FOUND - System has pugixml
#  PUGIXML_INCLUDES - The pugixml include directory (pugixml is header-only)
#  PUGIXML__LIBRARIES - The libraries needed to use pugixml
#
#  Set PUGIXML_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(PUGIXML_ROOT "/path/to/custom/pstreams") # prefer over system
#  find_package(pugixml)
#  ifPUGIXML_FOUND)
#    target_link_libraries (TARGET ${PUGIXML_LIBRARIES})
#  endif()

if(PUGIXML_INCLUDES)
  # Already in cache, be silent
  set (PUGIXML_FIND_QUIETLY TRUE)
endif()

FIND_PATH(PUGIXML_INCLUDES NAMES pugixml.hpp HINTS ${PUGIXML_ROOT}/include)
FIND_LIBRARY(PUGIXML_LIBRARIES NAMES pugixml HINTS ${PUGIXML_ROOT}/lib)

# Handle the QUIETLY and REQUIRED arguments and set PUGIXML_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(pugixml DEFAULT_MSG PUGIXML_LIBRARIES PUGIXML_INCLUDES)

MARK_AS_ADVANCED(PUGIXML_INCLUDES PUGIXML_LIBRARIES)
