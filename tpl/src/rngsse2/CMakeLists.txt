# CMake build for RNGSSE2
#
# Author: J. Bakosi
#
# Note: The library only contains C source, which compiles fine with clang,
# intel, and gnu compilers. However, pgicc complains about too long asm lines,
# which the pgc++ compiler builds fine, so we override cmake's guess, and build
# the c sources with the c++ compiler.

cmake_minimum_required(VERSION 2.8.5)

project(RNGSSELIB2 CXX)

include(GNUInstallDirs)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'RELEASE'")
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Build type." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

file(GLOB SOURCES "source/*.c")
file(GLOB HEADERS "include/*.h")

add_library(rngsse ${SOURCES})

# Override cmake's guess and build c sources with the c++ compiler
SET_SOURCE_FILES_PROPERTIES(${SOURCES} PROPERTIES LANGUAGE CXX)

INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
INSTALL(TARGETS rngsse
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)
