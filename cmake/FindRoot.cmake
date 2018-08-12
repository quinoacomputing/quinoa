################################################################################
#
# \file      cmake/FindRoot.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the Root library from CERN
#
################################################################################

# Find the Root library from CERN

#  Set ROOT_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(ROOT_ROOT "/path/to/custom/root") # prefer over system

#  find_package(ROOT COMPONENTS RIO Core Tree Hist)
#  if(ROOT_FOUND)
#    target_link_libraries(TARGET ${ROOT_LIBRARIES})
#  endif()


# If already in cache, be silent
if(ROOT_INCLUDE_DIRS AND ROOT_LIBRARIES)
  set (ROOT_FIND_QUIETLY TRUE)
endif()

### Find ROOT include files required
find_path(ROOT_INCLUDE_DIR NAMES TFile.h HINTS ${ROOT_ROOT}/include
                                               $ENV{ROOT_ROOT}/include
                                         PATH_SUFFIXES root root6)

### Find ROOT libraries required
set(ROOT_REQLIBS RIO Core Tree Hist Thread Net Imt Matrix MathCore)
foreach(lib ${ROOT_REQLIBS})
  find_library(ROOT_${lib}_LIBRARY NAMES ${lib}
               HINTS ${ROOT_ROOT}/lib $ENV{ROOT_ROOT}/lib
               PATH_SUFFIXES root root6)
  list(APPEND ROOT_LIBRARY "${ROOT_${lib}_LIBRARY}")
  list(APPEND ROOT_LIBRARY_VARS ROOT_${lib}_LIBRARY)
endforeach()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Root DEFAULT_MSG ${ROOT_LIBRARY_VARS} ROOT_INCLUDE_DIR)

if(NOT Root_FOUND)
  set(ROOT_INCLUDE_DIRS "")
  set(ROOT_LIBRARIES "")
else()
  set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})
  set(ROOT_LIBRARIES ${ROOT_LIBRARY})
endif()

MARK_AS_ADVANCED(ROOT_INCLUDE_DIRS ROOT_LIBRARIES)
