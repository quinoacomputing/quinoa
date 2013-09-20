# Build quinoa

cmake_minimum_required(VERSION 2.8.5)

project(QUINOA C CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Set Quinoa version number: major.minor
set(QUINOA_MAJOR_VER 0)
set(QUINOA_MINOR_VER 1)

# Set Quinoa release number
set(QUINOA_RELEASE "LA-CC-XX-XXX")

# Optionally turn off MKL calls
set(MKL_CALLS ON CACHE BOOL "Set to ON to call MKL functions")

IF(MKL_CALLS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_CALLS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_CALLS")
ENDIF(MKL_CALLS)

# Add compiler options needed for MKL
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")

# Set build type
set(QUINOA_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Set default third-party libraries directory
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
set(TPL_DIR ${PROJECT_SOURCE_DIR}/../tpl/install/${COMP} CACHE STRING
    "Third-party libraries directory")
MESSAGE(STATUS "Third-party library dir: " ${TPL_DIR})

# Find MPI and get its underlying compilers
include(GetMPICompilers)
get_mpi_compilers()

# Find OpenMP and set its compiler flag
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set compiler
set(QUINOA_MPI_COMPILER "${MPI_CXX_COMPILER}")
set(QUINOA_COMPILER "${UNDERLYING_CXX_COMPILER}")

# Grab current date to stick into the binary
execute_process(
  COMMAND "date" #"+%m/%d/%Y"
  OUTPUT_VARIABLE QUINOA_BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" QUINOA_BUILD_DATE "${QUINOA_BUILD_DATE}")

set(QUINOA_EXECUTABLE "quinoa")
set(QUINOA_RNGTEST_EXECUTABLE "rngtest")

# Get git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Add compiler flags
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-std=c++11 HAS_CPP0X)
if(HAS_CPP0X)
  add_definitions("-std=c++11")
  set (CMAKE_EXE_LINKER_FLAGS "-lstdc++")
endif()

check_cxx_compiler_flag(-Wall HAS_WALL)
if(HAS_WALL)
  add_definitions("-Wall")
endif()

check_cxx_compiler_flag(-pedantic HAS_PEDANTIC)
if(HAS_PEDANTIC)
  add_definitions("-pedantic")
  # suppress warning: ISO C++ forbids variable length array
  add_definitions("-Wno-vla")
endif()

check_cxx_compiler_flag(-Wremarks HAS_WREMARKS)
if(HAS_WREMARKS)
  add_definitions("-Wremarks")
endif()

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
  add_definitions("-Wno-vla-extension")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # Intel-specific compiler arguments:
  # Suppress the following Intel compiler remarks:
  #   111: statement is unreachable
  #   193: zero used for undefined preprocessing identifier
  #   304: access control not specified ("public" by default)
  #   383: value copied to temporary, reference to temporary used
  #   444: destructor for base class is not virtual
  #   981: operands are evaluated in unspecified order
  #  1418: external function definition with no prior declaration
  #  1599: declaration hides parameter
  #  3346: dynamic exception specifications are deprecated
  add_definitions("-diag-disable 111,193,304,383,444,981,1418,1599,3346")
  add_definitions("-opt-report=0")
endif()

# Make sure that the config file is in the search path
# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.h.in"
  "${PROJECT_BINARY_DIR}/Base/QuinoaConfig.h"
)

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

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/LinearAlgebra
                    ${QUINOA_SOURCE_DIR}/Main
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model/Mass/Beta
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/SLM
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/GLM
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model/Mix/Dirichlet
                    ${QUINOA_SOURCE_DIR}/Model/Mix/GenDirichlet
                    ${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Physics/HomMix
                    ${QUINOA_SOURCE_DIR}/Physics/HomHydro
                    ${QUINOA_SOURCE_DIR}/Physics/HomRT
                    ${QUINOA_SOURCE_DIR}/Physics/SPINSFlow
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Paradigm
                    ${TPL_DIR}/include
                    ${TPL_DIR}/include/pegtl
                    ${PROJECT_BINARY_DIR}/Base
                    ${MPI_INCLUDE_PATH}
)

# Components
add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(Geometry)
add_subdirectory(IO)
add_subdirectory(LinearAlgebra)
add_subdirectory(Mesh)
add_subdirectory(Model)
add_subdirectory(Paradigm)
add_subdirectory(Physics)
add_subdirectory(Random)
add_subdirectory(Statistics)
#add_subdirectory(FEM)
#add_subdirectory(Particles)
#add_subdirectory(SDE)
#add_subdirectory(WallTreatment)
#add_subdirectory(ThirdParty/exodiff)
#add_subdirectory(Utilities/forcenorm)
#add_subdirectory(Utilities/divset)
#add_subdirectory(test/rng)
add_subdirectory(Main)
