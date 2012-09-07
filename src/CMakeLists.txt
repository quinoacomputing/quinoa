cmake_minimum_required(VERSION 2.8.5)

project(QUINOA CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)

# Set Quinoa version number: major.minor
set(QUINOA_MAJOR_VERSION 0)
set(QUINOA_MINOR_VERSION 1)

# Set Quinoa release number
set(QUINOA_RELEASE "LA-CC-XX-XXX")

# Set platform string (from build script)
set(QUINOA_PLATFORM ${PLATFORM})

# Set third party prefix (from build script)
set(QUINOA_THIRD_PARTY_PREFIX ${THIRD_PARTY_PREFIX})

# Set build type
set(QUINOA_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Find OpenMP and set its compiler flag
#find_package(OpenMP REQUIRED)
#if(OPENMP_FOUND)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# Grab info about the compiler to stick into the binary
execute_process(
  COMMAND           ${CMAKE_CXX_COMPILER} "-dumpversion"
  OUTPUT_VARIABLE   CXX_COMPILER_VERSION
)
# Remove the newline character at the end
string(REGEX REPLACE "[\r\n]" "" CXX_COMPILER_VERSION
       "v${CXX_COMPILER_VERSION}")
# Set compiler string
set(QUINOA_COMPILER "${CMAKE_CXX_COMPILER} - ${CXX_COMPILER_VERSION}")

# Grab current date to stick into the binary
execute_process(
  COMMAND "date" "+%m/%d/%Y"
  OUTPUT_VARIABLE QUINOA_BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" QUINOA_BUILD_DATE "${QUINOA_BUILD_DATE}")

set (QUINOA_EXECUTABLE "quinoa")

# Set DEBUG bool depending on CMAKE_BUILD_TYPE
IF(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set (DEBUG OFF)
ELSE(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set (DEBUG ON)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "RELEASE")

# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.h.in"
  "${PROJECT_BINARY_DIR}/Base/QuinoaConfig.h"
)

# Add compiler flags
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wall HAS_WALL)
IF(HAS_WALL)
  add_definitions("-Wall")
ENDIF(HAS_WALL)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++0x HAS_CPP0x)
IF(HAS_CPP0x)
  add_definitions("-std=c++0x")
ENDIF(HAS_CPP0x)

### Set RPATH so that shared libs paths stay in the binary after install ###
# Don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# When building, don't use the install RPATH
# (use later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")


# Make sure that the config file is in the search path
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/Base")

# Infrastructure
add_subdirectory(Base)
#add_subdirectory(Control)
#add_subdirectory(FEM)
add_subdirectory(IO)
add_subdirectory(LinearAlgebra)
add_subdirectory(Mesh)
#add_subdirectory(ParseControl)
#add_subdirectory(Particles)
#add_subdirectory(SDE)
#add_subdirectory(WallTreatment)
#add_subdirectory(ThirdParty/exodiff)
#add_subdirectory(Utilities/forcenorm)
#add_subdirectory(Utilities/divset)

# Main
add_subdirectory(Main)
