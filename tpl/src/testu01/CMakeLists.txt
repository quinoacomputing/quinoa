# CMake build for TestU01
#
# Author: J. Bakosi
#
# Minimum cmake file for building TestU01.
#
# Notes:
#  * Of the original TestU01 distribution we only build the testu01, probdist,
#    and mylib libraries.
#  * The header files for TestU01 are generated from .tex files during running
#    make in the original TestU01 distribution preceded by configure. We simply
#    take over those header files generated and build the libraries with cmake.
#  * Rationale: autotools crapped out with cmake's external_project_add from a
#    SOURCE_DIR specified (only worked from a tarball)

cmake_minimum_required(VERSION 2.8.5)

project(TESTU01 C)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'DEBUG'")
  set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Add compiler flags
include(CheckCCompilerFlag)

macro (require_compiler_flag flag)
  check_c_compiler_flag(${flag} COMPILER_HAS${flag})
  if(COMPILER_HAS${flag})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
  else()
    message(FATAL "The compiler must support the -Wno-return-type argument")
  endif()
endmacro (require_compiler_flag)

# Require the following arguments for all compilers
require_compiler_flag("-Wno-return-type")

file(GLOB PROBDIST_SOURCES "probdist/*.c")
file(GLOB MYLIB_SOURCES "mylib/*.c")
file(GLOB TESTU01_SOURCES "testu01/*.c")
file(GLOB HEADERS "include/*.h")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(probdist ${PROBDIST_SOURCES})
add_library(mylib ${MYLIB_SOURCES})
add_library(testu01 ${TESTU01_SOURCES})

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(TARGETS probdist
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development)
INSTALL(TARGETS mylib
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development)
INSTALL(TARGETS testu01
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development)
