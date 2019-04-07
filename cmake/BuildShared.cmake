################################################################################
#
# \file      BuildShared.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Set default value for building shared libs if none was specified
#
################################################################################

if(NOT DEFINED BUILD_SHARED_LIBS)

  # Detect a Cray machine. This is based on testing for the following
  # environment variables. At least one of these is always defined on a Cray,
  # depending on what programming environment is loaded among the modules.
  if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
      NOT DEFINED ENV{CRAY_PRGENVGNU} AND
      NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
      NOT DEFINED ENV{CRAY_PRGENVINTEL})
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries. Possible values: ON | OFF")
    message(STATUS "BUILD_SHARED_LIBS not specified, setting to 'ON'")
  else()
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries. Possible values: ON | OFF")
    message(STATUS "Cray detected, setting BUILD_SHARED_LIBS to 'OFF'")
  endif()

endif()

message(STATUS "BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
