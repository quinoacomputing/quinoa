################################################################################
#
# \file      cmake/FindPEGTL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find PEGTL
#
################################################################################

# See PEGTL: https://github.com/ColinH/PEGTL
#
#  PEGTL_FOUND - System has PEGTL
#  PEGTL_INCLUDE_DIRS - The PEGTL include directory
#
#  Set the PEGTL_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(PEGTL_ROOT "/path/to/custom/pegtl") # prefer over system
#  find_package(PEGTL)
#  include_directories(${PEGTL_INCLUDE_DIRS})

# If already in cache, be silent
if(PEGTL_INCLUDE_DIRS)
  set (PEGTL_FIND_QUIETLY TRUE)
endif()

find_path(PEGTL_INCLUDE_DIR NAMES pegtl.hpp
                            HINTS ${PEGTL_ROOT}/include
                                  $ENV{PEGTL_ROOT}
                            PATH_SUFFIXES pegtl tao/pegtl pegtl/include/tao)

set(PEGTL_INCLUDE_DIRS ${PEGTL_INCLUDE_DIR})

# Handle the QUIETLY and REQUIRED arguments and set PEGTL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PEGTL DEFAULT_MSG PEGTL_INCLUDE_DIRS)

MARK_AS_ADVANCED(PEGTL_INCLUDE_DIRS)
