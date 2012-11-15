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

# Set link strategy (from build script)
set(LINK_STRATEGY ${LINK_STRATEGY}
    CACHE STRING "DYNAMIC, STATIC (Set this at the command line)"
    FORCE)

# Add compiler options needed for MKL
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")

# Make sure ICC uses gcc-4.7 compatibililty
if(PLATFORM MATCHES "linux-intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -gcc-version=470 -gcc-name=gcc-4.7")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gcc-version=470 -gxx-name=g++-4.7")
endif()

# Make sure ICC finds the gcc-4.7 libs, installed in custom location on darwin
if(PLATFORM MATCHES "darwin-intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -gcc-version=470 -gcc-name=gcc -cxxlib=/home/opt/gcc/4.7.1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gcc-version=470 -gxx-name=g++ -cxxlib=/home/opt/gcc/4.7.1")
endif()

# Set third party prefix (from build script)
set(QUINOA_THIRD_PARTY_PREFIX ${THIRD_PARTY_PREFIX})

# Set build type
set(QUINOA_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Find OpenMP and set its compiler flag
find_package(OpenMP)
if(OPENMP_FOUND AND NOT PLATFORM MATCHES "clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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
  COMMAND "date" #"+%m/%d/%Y"
  OUTPUT_VARIABLE QUINOA_BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" QUINOA_BUILD_DATE "${QUINOA_BUILD_DATE}")

set (QUINOA_EXECUTABLE "quinoa")

# Set DEBUG bool depending on CMAKE_BUILD_TYPE
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set (DEBUG off)
else()
  set (DEBUG on)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModule/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Make sure that the config file is in the search path
# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.h.in"
  "${PROJECT_BINARY_DIR}/Base/QuinoaConfig.h"
)

# Add compiler flags
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-std=c++0x HAS_CPP0X)
if(HAS_CPP0X)
  add_definitions("-std=c++0x")
  set (CMAKE_EXE_LINKER_FLAGS "-lstdc++")
endif()

check_cxx_compiler_flag(-Wall HAS_WALL)
if(HAS_WALL)
  add_definitions("-Wall")
endif()

check_cxx_compiler_flag(-Wremarks HAS_WREMARKS)
if(HAS_WREMARKS)
  # Turn on remarks
  add_definitions("-Wremarks")
  # Disregard the following remarks:
  #   981: operands are evaluated in unspecified order
  #   383: value copied to temporary, reference to temporary used
  add_definitions("-diag-disable 981")
  add_definitions("-diag-disable 383")
endif()

check_cxx_compiler_flag(-opt-report HAS_OPTREPORT)
if(HAS_OPTREPORT)
  add_definitions("-opt-report=0")
endif()

### Set RPATH so that shared libs paths stay in the binary after install ###
# Don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH false)

# When building, don't use the install RPATH
# (use later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH false)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/Base")

# Infrastructure
add_subdirectory(Base)
add_subdirectory(IO)
add_subdirectory(LinearAlgebra)
add_subdirectory(Mesh)
add_subdirectory(Paradigm)
add_subdirectory(Random)
add_subdirectory(Statistics)
add_subdirectory(Model)
add_subdirectory(Model/MixModel)
add_subdirectory(Model/HomDirichlet)
add_subdirectory(Model/HomGenDirichlet)
#add_subdirectory(FEM)
#add_subdirectory(Particles)
#add_subdirectory(SDE)
#add_subdirectory(WallTreatment)
#add_subdirectory(ThirdParty/exodiff)
#add_subdirectory(Utilities/forcenorm)
#add_subdirectory(Utilities/divset)

# Main
add_subdirectory(Main)
