################################################################################
#
# \file      cmake/FindRoot.cmake
# \author    A. Pakki
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Root library from CERN
#
################################################################################

# Find the Root library from CERN
#
# If already in cache, be silent
if(ROOT_INCLUDE_DIRS AND ROOT_LIBRARIES)
  set (ROOT_FIND_QUIETLY TRUE)
endif()

find_path(ROOT_INCLUDE_DIR NAMES TFile.h)

### Find ROOT libraries required
find_library(ROOT_RIO_LIBRARY NAMES RIO)
find_library(ROOT_CORE_LIBRARY NAMES Core)
find_library(ROOT_TREE_LIBRARY NAMES Tree)
find_library(ROOT_GRAPH_LIBRARY NAMES Hist)
find_library(ROOT_THREAD_LIBRARY NAMES Thread)
find_library(ROOT_NET_LIBRARY NAMES Net)
find_library(ROOT_IMT_LIBRARY NAMES Imt)
find_library(ROOT_MATRIX_LIBRARY NAMES Matrix)
find_library(ROOT_MATHCORE_LIBRARY NAMES MathCore)

### Link the libraries as one
set(ROOT_LIBRARY ${ROOT_RIO_LIBRARY}
                 ${ROOT_CORE_LIBRARY}
                 ${ROOT_TREE_LIBRARY}
                 ${ROOT_GRAPH_LIBRARY}
                 ${ROOT_THREAD_LIBRARY}
                 ${ROOT_NET_LIBRARY}
                 ${ROOT_IMT_LIBRARY}
                 ${ROOT_MATRIX_LIBRARY}
                 ${ROOT_MATHCORE_LIBRARY})

set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})
set(ROOT_LIBRARIES ${ROOT_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Root DEFAULT_MSG ROOT_LIBRARIES ROOT_INCLUDE_DIRS)

if(NOT Root_FOUND)
  set(ROOT_INCLUDE_DIRS "")
  set(ROOT_LIBRARIES "")
endif()

MARK_AS_ADVANCED(ROOT_INCLUDE_DIRS ROOT_LIBRARIES)
