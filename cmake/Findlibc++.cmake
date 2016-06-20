################################################################################
#
# \file      cmake/Findlibc++.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find libc++
# \date      Mon 20 Jun 2016 12:58:37 PM MDT
#
################################################################################

# Find libc++.
# See libc++: http://libcxx.llvm.org, libc++abi: http://libcxxabi.llvm.org.
#
#  libc++_FOUND - System has libc++
#  libc++_INCLUDES - The libc++ include directory
#  libc++_LIBRARIES - The libraries needed to use libc++
#  libc++abi_LIBRARIES - The libraries needed to use libc++abi
#
#  Set the LIBCXX_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(LIBCXX_ROOT "/path/to/custom/libc++") # prefer over system
#  find_package(libc++)
#  if(libc++_FOUND)
#    target_link_libraries (TARGET ${libc++_LIBRARIES} ${libc++abi_LIBRARIES})
#  endif()

# If already in cache, be silent
if(libc++_INCLUDES AND libc++_LIBRARIES AND libc++abi_LIBRARIES)
  set (libc++_FIND_QUIETLY TRUE)
endif()

find_path(libc++_INCLUDES NAMES cxxabi.h HINTS ${LIBCXX_ROOT}/include
                                               /usr/include/c++/v1
                                               $ENV{LIBCXX_ROOT}/include/c++/v1)

if(BUILD_SHARED_LIBS)
  find_library(libc++_LIBRARIES NAMES c++ HINTS ${LIBCXX_ROOT}/lib
                                                $ENV{LIBCXX_ROOT}/lib)
  find_library(libc++abi_LIBRARIES NAMES c++abi HINTS ${LIBCXX_ROOT}/lib
                                                $ENV{LIBCXX_ROOT}/lib)
else()
  find_library(libc++_LIBRARIES NAMES libc++.a HINTS ${LIBCXX_ROOT}/lib
                                                     $ENV{LIBCXX_ROOT}/lib)
  find_library(libc++abi_LIBRARIES NAMES libc++abi.a HINTS ${LIBCXX_ROOT}/lib
                                                     $ENV{LIBCXX_ROOT}/lib)
endif()

# Handle the QUIETLY and REQUIRED arguments and set libc++_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(libc++ DEFAULT_MSG libc++_LIBRARIES libc++abi_LIBRARIES libc++_INCLUDES)

MARK_AS_ADVANCED(libc++_INCLUDES libc++_LIBRARIES libc++abi_LIBRARIES)
