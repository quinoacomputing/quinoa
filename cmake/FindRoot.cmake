################################################################################
#
# \file      cmake/FindHypre.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Root library from CERN
# \date      Fri 20 Jan 2017 12:28:14 PM MST
#
################################################################################

# Find the Root library from CERN
#
# If already in cache, be silent
if(ROOT_INCLUDE_DIRS AND ROOT_LIBRARIES)
  set (ROOT_FIND_QUIETLY TRUE)
endif()

find_path(ROOT_INCLUDE_DIR NAMES TFile.h 
                           HINTS ${ROOTSYS}/include
)

find_library(ROOT_LIBRARY NAMES RIO Core TTree
	                  HINTS ${ROOTSYS}/lib 
)

set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})
set(ROOT_LIBRARIES ${ROOT_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Root DEFAULT_MSG ROOT_LIBRARIES ROOT_INCLUDE_DIRS)

MARK_AS_ADVANCED(ROOT_INCLUDE_DIRS ROOT_LIBRARIES)
