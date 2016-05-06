################################################################################
#
# \file      cmake/BuildShared.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Set default value for building shared libs if none was specified
# \date      Fri 06 May 2016 06:40:37 AM MDT
#
################################################################################

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries. Possible values: ON | OFF")
if(NOT DEFINED BUILD_SHARED_LIBS)
  message(STATUS "BUILD_SHARED_LIBS not specified, setting to 'ON'")
else()
  if(NOT BUILD_SHARED_LIBS AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Static linking only supported with the GNU compiler suite, overriding BUILD_SHARED_LIBS = off -> on")
    set(BUILD_SHARED_LIBS on)
  endif()
  message(STATUS "BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
endif()
