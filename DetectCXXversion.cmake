################################################################################
#
# \file      cmake/DetectCXXversion.cmake
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Detect C++ compiler major, minor, and patch version
#
################################################################################

string(REGEX MATCH "([0-9]*)\\.([0-9]*)\\.([0-9]*)" major ${CMAKE_CXX_COMPILER_VERSION})
set(CMAKE_CXX_COMPILER_MAJOR ${CMAKE_MATCH_1})
set(CMAKE_CXX_COMPILER_MINOR ${CMAKE_MATCH_2})
set(CMAKE_CXX_COMPILER_PATCH ${CMAKE_MATCH_3})
