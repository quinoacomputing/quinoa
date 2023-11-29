################################################################################
#
# \file      FindLIBCXX.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Find libc++
#
################################################################################

# libc++: http://libcxx.llvm.org, libc++abi: http://libcxxabi.llvm.org
#
#  LIBCXX_FOUND - System has libc++
#  LIBCXX_INCLUDE_DIRS - The libc++ include directory
#  LIBCXX_LIBRARIES - The libraries needed to use libc++
#  LIBCXXABI_LIBRARIES - The libraries needed to use libc++abi
#
#  Set the LIBCXX_PATH cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(LIBCXX_PATH "/path/to/custom/libc++") # prefer over system
#  find_package(LIBCXX)
#  if(LIBCXX_FOUND)
#    target_link_libraries (TARGET ${LIBCXX_LIBRARIES} ${LIBCXXABI_LIBRARIES})
#  endif()

# If already in cache, be silent
if(LIBCXX_INCLUDE_DIRS AND LIBCXX_LIBRARIES AND LIBCXXABI_LIBRARIES)
  set (LIBCXX_FIND_QUIETLY TRUE)
endif()

find_path(LIBCXX_INCLUDE_DIR NAMES cxxabi.h HINTS ${LIBCXX_PATH}/include
                                            /usr/include/c++/v1
                                            $ENV{LIBCXX_PATH}/include/c++/v1)

find_library(LIBCXX_LIBRARY NAMES c++ HINTS ${LIBCXX_PATH}/lib
                                            $ENV{LIBCXX_PATH}/lib)
find_library(LIBCXXABI_LIBRARY NAMES c++abi HINTS ${LIBCXX_PATH}/lib
                                            $ENV{LIBCXX_PATH}/lib)

if(LIBCXX_INCLUDE_DIR AND LIBCXX_LIBRARY AND LIBCXXABI_LIBRARY)
  set(LIBCXX_INCLUDE_DIRS ${LIBCXX_INCLUDE_DIR})
  set(LIBCXX_LIBRARIES ${LIBCXX_LIBRARY})
  set(LIBCXXABI_LIBRARIES ${LIBCXXABI_LIBRARY})
else()
  set(LIBCXX_INCLUDE_DIRS "")
  set(LIBCXX_LIBRARIES "")
  set(LIBCXXABI_LIBRARIES "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set LIBCXX_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBCXX DEFAULT_MSG LIBCXX_LIBRARIES LIBCXXABI_LIBRARIES LIBCXX_INCLUDE_DIRS)

MARK_AS_ADVANCED(LIBCXX_INCLUDE_DIRS LIBCXX_LIBRARIES LIBCXXABI_LIBRARIES)
