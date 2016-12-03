#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find H5Part
# Find the native H5Part headers and libraries.
#
#  H5Part_INCLUDE_DIRS - where to find flecsi.h, etc.
#  H5Part_LIBRARIES    - List of libraries when using flecsi.
#  H5Part_FOUND        - True if flecsi found.

# Look for the header file.
FIND_PATH(H5Part_INCLUDE_DIR NAMES H5Part.h)

# Look for the library.
if(NOT BUILD_SHARED_LIBS)
  FIND_LIBRARY(H5Part_LIBRARY NAMES libH5Part.a)
else()
  FIND_LIBRARY(H5Part_LIBRARY NAMES H5Part)
endif()

# handle the QUIETLY and REQUIRED arguments and set H5Part_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(H5Part H5Part_LIBRARY H5Part_INCLUDE_DIR)

# Copy the results to the output variables.
SET(H5Part_LIBRARIES ${H5Part_LIBRARY})
SET(H5Part_INCLUDE_DIRS ${H5Part_INCLUDE_DIR})

MARK_AS_ADVANCED(H5Part_INCLUDE_DIR H5Part_LIBRARY)
