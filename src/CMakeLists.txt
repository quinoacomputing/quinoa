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

# Get MPI & its underlying compilers
find_package(MPI REQUIRED)

# Find out underlying C compiler
execute_process(
  COMMAND           ${MPI_C_COMPILER} "-showme:command"
  OUTPUT_VARIABLE   UNDERLYING_C_COMPILER 
)
# Remove the newline character at the end of UNDERLYING_C_COMPILER
string(REGEX REPLACE "[\r\n]" "" UNDERLYING_C_COMPILER "${UNDERLYING_C_COMPILER}")

# Find out underlying C++ compiler
execute_process(
  COMMAND           ${MPI_CXX_COMPILER} "-showme:command"
  OUTPUT_VARIABLE   UNDERLYING_CXX_COMPILER 
)
# Remove the newline character at the end of UNDERLYING_CXX_COMPILER
string(REGEX REPLACE "[\r\n]" "" UNDERLYING_CXX_COMPILER "${UNDERLYING_CXX_COMPILER}")

# Echo compilers
MESSAGE(STATUS "MPI C compiler: " ${MPI_C_COMPILER})
MESSAGE(STATUS "MPI C++ compiler: " ${MPI_CXX_COMPILER})
MESSAGE(STATUS "MPI underlying C compiler: " ${UNDERLYING_C_COMPILER})
MESSAGE(STATUS "MPI underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})

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

# Mass model options
# Available mass models
set(MASS_VALUES "BETA")
# Initialize all to off
set(QUINOA_BETA off)  # 0
# Set default and select from list
set(QUINOA_MODEL_MASS "BETA"
    CACHE STRING "Select mass model. Default: BETA (Beta distribution model). Available options: ${MASS_VALUES}")
SET_PROPERTY (CACHE QUINOA_MODEL_MASS PROPERTY STRINGS ${MASS_VALUES})
STRING (TOUPPER ${QUINOA_MODEL_MASS} QUINOA_MODEL_MASS)
LIST (FIND MASS_VALUES ${QUINOA_MODEL_MASS} MASS_INDEX)
# Evaluate selected model and put in define for it
IF (${MASS_INDEX} EQUAL 0)
  set(QUINOA_BETA on)
ELSEIF (${MASS_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Mass model ${QUINOA_MODEL_MASS} not supported, valid entries are ${MASS_VALUES}")
ENDIF()

# Hydro model options
# Available hydrodyanmics models
set(HYDRO_VALUES "SLM" "GLM")
# Initialize all to off
set(QUINOA_SLM off)  # 0
set(QUINOA_GLM off)  # 1
# Set default and select from list
set(QUINOA_MODEL_HYDRO "SLM"
    CACHE STRING "Select hydrodynamics model. Default: SLM (Simplified Langevin model). Available options: ${HYDRO_VALUES}")
SET_PROPERTY (CACHE QUINOA_MODEL_HYDRO PROPERTY STRINGS ${HYDRO_VALUES})
STRING (TOUPPER ${QUINOA_MODEL_HYDRO} QUINOA_MODEL_HYDRO)
LIST (FIND HYDRO_VALUES ${QUINOA_MODEL_HYDRO} HYDRO_INDEX)
# Evaluate selected model and put in define for it
IF (${HYDRO_INDEX} EQUAL 0)
  set(QUINOA_SLM on)
ELSEIF (${HYDRO_INDEX} EQUAL 1)
  set(QUINOA_GLM on)
ELSEIF (${HYDRO_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Hydrodynamics model ${QUINOA_MODEL_HYDRO} not supported, valid entries are ${HYDRO_VALUES}")
ENDIF()

# Mix model options
# Available mix models
set(MIX_VALUES "DIR" "GENDIR")
# Initialize all to off
set(QUINOA_DIRICHLET off)              # 0
set(QUINOA_GENERALIZED_DIRICHLET off)  # 1
# Set default and select from list
set(QUINOA_MODEL_MIX "DIR"
    CACHE STRING "Select mixing model. Default: DIR (Dirichlet model). Available options: ${MIX_VALUES}")
SET_PROPERTY (CACHE QUINOA_MODEL_MIX PROPERTY STRINGS ${MIX_VALUES})
STRING (TOUPPER ${QUINOA_MODEL_MIX} QUINOA_MODEL_MIX)
LIST (FIND MIX_VALUES ${QUINOA_MODEL_MIX} MIX_INDEX)
# Evaluate selected model and put in define for it
IF (${MIX_INDEX} EQUAL 0)
  set(QUINOA_DIRICHLET on)
ELSEIF (${MIX_INDEX} EQUAL 1)
  set(QUINOA_GENERALIZED_DIRICHLET on)
ELSEIF (${MIX_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Mix model ${QUINOA_MODEL_MIX} not supported, valid entries are ${MIX_VALUES}")
ENDIF()

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
  #   304: access control not specified ("public" by default)
  #   383: value copied to temporary, reference to temporary used
  #   981: operands are evaluated in unspecified order
  #  1418: external function definition with no prior declaration
  #  1599: declaration hides parameter
  add_definitions("-diag-disable 111,304,383,981,1418,1599")
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

include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/Base)
include_directories(${MPI_INCLUDE_PATH})

# Components
add_subdirectory(Base)
#add_subdirectory(Geometry)
add_subdirectory(IO)
add_subdirectory(LinearAlgebra)
add_subdirectory(Mesh)
add_subdirectory(Model/Mass/Beta)
add_subdirectory(Model/Mix/Dirichlet)
add_subdirectory(Model/Mix/GeneralizedDirichlet)
add_subdirectory(Model/Hydro/SimplifiedLangevin)
add_subdirectory(Model/Hydro/GeneralizedLangevin)
add_subdirectory(Paradigm)
add_subdirectory(Parser)
add_subdirectory(Physics)
add_subdirectory(Physics/HomMix)
add_subdirectory(Physics/HomHydro)
add_subdirectory(Physics/HomRT)
add_subdirectory(Physics/SPINSFlow)
add_subdirectory(Random)
add_subdirectory(Statistics)
#add_subdirectory(FEM)
#add_subdirectory(Particles)
#add_subdirectory(SDE)
#add_subdirectory(WallTreatment)
#add_subdirectory(ThirdParty/exodiff)
#add_subdirectory(Utilities/forcenorm)
#add_subdirectory(Utilities/divset)
add_subdirectory(Main)
