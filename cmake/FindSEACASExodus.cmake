#------------------------------------------------------------------------------#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#------------------------------------------------------------------------------#

# - Find libexoIIv2c
# Find the native SEACASExodus headers and libraries.
#
#  SEACASExodus_INCLUDE_DIRS - where to find exodusII.h, etc.
#  SEACASExodus_LIBRARIES    - List of libraries when using xoIIv2c.
#  SEACASExodus_FOUND        - True if exodus found.
#

find_path(SEACASExodus_INCLUDE_DIR exodusII.h PATH_SUFFIXES trilinos
          HINTS $ENV{EXODUS_ROOT}/include)

#v6.09 calls it libexodus
#debian calls it libexoIIv2
#other distros libexoIIv2c
find_library(SEACASExodus_LIBRARY NAMES SEACASExodus exodus exoIIv2 exoIIv2c
                                  HINTS $ENV{EXODUS_ROOT}/lib
                                  PATH_SUFFIXES trilinos)

set(SEACASExodus_LIBRARIES ${SEACASExodus_LIBRARY})
set(SEACASExodus_INCLUDE_DIRS ${SEACASExodus_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SEACASExodus_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SEACASExodus DEFAULT_MSG SEACASExodus_LIBRARY SEACASExodus_INCLUDE_DIR )

mark_as_advanced(SEACASExodus_INCLUDE_DIR SEACASExodus_LIBRARY)
