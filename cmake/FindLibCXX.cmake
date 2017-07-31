################################################################################
#
# \file      cmake/FindLibCXX.cmake
# \copyright 2016-2017, Los Alamos National Security, LLC.
# \brief     Find libc++
#
################################################################################

# Find libc++.
# See libc++: http://libcxx.llvm.org, libc++abi: http://libcxxabi.llvm.org.
#
#  LIBCXX_FOUND - System has libc++
#  LIBCXX_INCLUDE_DIRS - The libc++ include directory
#  LIBCXX_LIBRARIES - The libraries needed to use libc++
#  LIBCXXABI_LIBRARIES - The libraries needed to use libc++abi
#
#  Set the LIBCXX_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(LIBCXX_ROOT "/path/to/custom/libc++") # prefer over system
#  find_package(LibCXX)
#  if(LIBCXX_FOUND)
#    target_link_libraries (TARGET ${LIBCXX_LIBRARIES} ${LIBCXXABI_LIBRARIES})
#  endif()

# If already in cache, be silent
if(LIBCXX_INCLUDE_DIRS AND LIBCXX_LIBRARIES AND LIBCXXABI_LIBRARIES)
  set (LIBCXX_FIND_QUIETLY TRUE)
endif()

find_path(LIBCXX_INCLUDE_DIR NAMES cxxabi.h HINTS ${LIBCXX_ROOT}/include
                                            /usr/include/c++/v1
                                            $ENV{LIBCXX_ROOT}/include/c++/v1)

if(BUILD_SHARED_LIBS)
  find_library(LIBCXX_LIBRARIES NAMES c++ HINTS ${LIBCXX_ROOT}/lib
                                                $ENV{LIBCXX_ROOT}/lib)
  find_library(LIBCXXABI_LIBRARIES NAMES c++abi HINTS ${LIBCXX_ROOT}/lib
                                                $ENV{LIBCXX_ROOT}/lib)
else()
  find_library(LIBCXX_LIBRARIES NAMES libc++.a HINTS ${LIBCXX_ROOT}/lib
                                                     $ENV{LIBCXX_ROOT}/lib)
  if(ARCH MATCHES "ppc64")
    set(LIBCXXABI_LIBRARIES " ")
  else()
    find_library(LIBCXXABI_LIBRARIES NAMES libc++abi.a HINTS ${LIBCXX_ROOT}/lib
                                                       $ENV{LIBCXX_ROOT}/lib)
  endif()
endif()

set(LIBCXX_INCLUDE_DIRS ${LIBCXX_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set LIBCXX_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBCXX DEFAULT_MSG LIBCXX_LIBRARIES LIBCXXABI_LIBRARIES LIBCXX_INCLUDE_DIRS)

MARK_AS_ADVANCED(LIBCXX_INCLUDE_DIRS LIBCXX_LIBRARIES LIBCXXABI_LIBRARIES)
