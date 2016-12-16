################################################################################
#
# \file      cmake/FindPEGTL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find PEGTL
# \date      Fri 16 Dec 2016 08:36:38 AM MST
#
################################################################################

# See PEGTL: https://github.com/ColinH/PEGTL
#
#  PEGTL_FOUND - System has PEGTL
#  PEGTL_INCLUDE_PATH - The PEGTL include directory
#
#  Set the PEGTL_ROOT cmake variable or shell environment variable before
#  calling find_package to a path to add an additional search path, e.g.,
#
#  Usage:
#
#  set(PEGTL_ROOT "/path/to/custom/pegtl") # prefer over system
#  find_package(PEGTL)
#  include_directories(${PEGTL_INCLUDE_PATH})

# If already in cache, be silent
if(PEGTL_INCLUDE_PATH)
  set (PEGTL_FIND_QUIETLY TRUE)
endif()

find_path(PEGTL_INCLUDE_PATH NAMES pegtl.hh
                             HINTS ${PEGTL_ROOT}/include
                             PATH_SUFFIXES pegtl)

# Handle the QUIETLY and REQUIRED arguments and set PEGTL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PEGTL DEFAULT_MSG PEGTL_INCLUDE_PATH)

MARK_AS_ADVANCED(PEGTL_INCLUDE_PATH)
