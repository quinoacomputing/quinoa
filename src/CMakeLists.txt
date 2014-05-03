# Build quinoa

cmake_minimum_required(VERSION 2.8.5)

project(QUINOA C CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
disallow_in_source_builds()

# Set version number: major.minor
set(MAJOR_VER 0)
set(MINOR_VER 1)

# Set release number
set(RELEASE "LA-CC-XX-XXX")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'DEBUG'")
  set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type. Possible values: DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Set build type
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Build shared-libs by default
set(BUILD_SHARED_LIBS on)

# If not specified on the command-line, set default third-party libraries
# directory based on cmake-supplied compiler ID
if (NOT TPL_DIR)
  string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
  set(REL_TPL_DIR ${PROJECT_SOURCE_DIR}/../tpl/install/${COMP})
  get_filename_component(ABS_TPL_DIR ${REL_TPL_DIR} ABSOLUTE)
  set(TPL_DIR ${ABS_TPL_DIR} CACHE STRING "Third-party libraries directory.")
endif()

# Echo status on third-party libraries directory
if (TPL_DIR)
  MESSAGE(STATUS "Third-party library dir: " ${TPL_DIR})
else()
  message(FATAL_ERROR "Third-party library dir '${TPL_DIR}' is not specified")
endif()

# Find MPI and get its underlying compilers
include(MPICompilers)
get_mpi_compilers()

# Set compilers
set(COMPILER ${UNDERLYING_CXX_COMPILER})
set(MPI_COMPILER ${MPI_CXX_COMPILER})
set(CHARM_COMPILER ${TPL_DIR}/charm/bin/charmc)
set(AMPI_COMPILER ${TPL_DIR}/charm/bin/ampicxx)
message(STATUS "Charm++ compiler: " ${CHARM_COMPILER})
message(STATUS "Charm++/AMPI compiler: " ${AMPI_COMPILER})

set(CMAKE_C_COMPILER ${CHARM_COMPILER} CACHE STRING "C compiler" FORCE)
set(CMAKE_CXX_COMPILER ${CHARM_COMPILER} CACHE STRING "CXX compiler" FORCE)

# Find OpenMP and set its compiler flag
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Query build hostname
site_name(HOSTNAME)

# Grab current date to stick into the binary
execute_process(
  COMMAND "date" #"+%m/%d/%Y"
  OUTPUT_VARIABLE BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

# Set names of executables
set(QUINOA_EXECUTABLE "quinoa")
set(RNGTEST_EXECUTABLE "rngtest")
set(MESHCONV_EXECUTABLE "meshconv")

# Get git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Add compiler flags
include(CheckCXXCompilerFlag)

macro (add_compiler_flag flag)
  check_cxx_compiler_flag(${flag} COMPILER_HAS${flag})
  if(COMPILER_HAS${flag})
    add_definitions(${flag})
  endif()
endmacro (add_compiler_flag)

# Try and use the following arguments for all compilers
add_compiler_flag("-Wall")
add_compiler_flag("-pedantic")
add_compiler_flag("-Wextra")
add_compiler_flag("-Wshadow")
add_compiler_flag("-Wformat")
add_compiler_flag("-Wsign-compare")
add_compiler_flag("-Wpointer-arith")
add_compiler_flag("-Winvalid-pch")
add_compiler_flag("-Wcast-align")
add_compiler_flag("-Wwrite-strings")
add_compiler_flag("-Wunused")
add_compiler_flag("-Woverloaded-virtual")
add_compiler_flag("-ansi")
add_compiler_flag("-fdiagnostics-show-option")

#add_compiler_flag("-default-to-aout")

# Suppress warning for all compilers: ISO C++ forbids variable length array
add_compiler_flag("-Wno-vla")

check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CPP11)
if(COMPILER_HAS_CPP11)
  add_definitions("-std=c++11")
else()
  MESSAGE(FATAL_ERROR "Need a compiler with C++11 support.")
endif()

# Include function that allows to set standard C++ library
include(libstdcxx)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")  # clang-specific settings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")

  # Set libc++ as the default standard C++ library and allow selection
  setlibstdcxx("libc++")

  # Experimental: use the C++14 standard
  add_definitions("-std=c++1y")

  # Don't warn on variable-length arrays
  add_compiler_flag("-Wno-vla-extension")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings
  # Set libstdc++ as the default standard C++ library and allow selection
  setlibstdcxx("libstdc++")

  add_compiler_flag("-Wno-literal-suffix")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings
  add_compiler_flag("-Wremarks")
  # Suppress the following Intel compiler remarks:
  #   111: statement is unreachable
  #   193: zero used for undefined preprocessing identifier
  #   304: access control not specified ("public" by default)
  #   383: value copied to temporary, reference to temporary used
  #   444: destructor for base class is not virtual
  #   869: parameter was never referenced
  #   981: operands are evaluated in unspecified order
  #  1418: external function definition with no prior declaration
  #  2415: variable of static storage duration was declared but never referenced
  #  3346: dynamic exception specifications are deprecated
  add_definitions("-diag-disable 111,193,304,383,444,869,981,1418,2415,3346")
  add_definitions("-opt-report=0")
  #add_compiler_flag("-openmp-report")

  # Link to thread-safe MPI library
  add_compiler_flag("-mt_mpi")

  # Set libstdc++ as the default standard C++ library and allow selection
  setlibstdcxx("libstdc++")

endif()

# Set compiler flags for use of libc++ (If nothing is given, libstdc++ is the
# default thus does not require anything special.)
if(LIBCXX_CPP)
   # Use libc++ (http://libcxx.llvm.org)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
   # Use libc++abi (http://libcxxabi.llvm.org)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lc++abi")
endif()

# Set data layout for particle properties
# Available options
set(LAYOUT_VALUES "particle" "equation")
# Initialize all to off
set(LAYOUT_PARTICLE_MAJOR off)  # 0
set(LAYOUT_EQUATION_MAJOR off)  # 1
# Set default and select from list
set(LAYOUT "particle" CACHE STRING "Particle properties layout. Default: (particle-major). Available options: ${LAYOUT_VALUES}(-major).")
SET_PROPERTY (CACHE LAYOUT PROPERTY STRINGS ${LAYOUT_VALUES})
STRING (TOLOWER ${LAYOUT} LAYOUT)
LIST (FIND LAYOUT_VALUES ${LAYOUT} LAYOUT_INDEX)
# Evaluate selected option and put in a define for it
IF (${LAYOUT_INDEX} EQUAL 0)
  set(LAYOUT_PARTICLE_MAJOR on)
ELSEIF (${LAYOUT_INDEX} EQUAL 1)
  set(LAYOUT_EQUATION_MAJOR on)
ELSEIF (${LAYOUT_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Particle properties layout '${LAYOUT}' not supported, valid entries are ${LAYOUT_VALUES}(-major).")
ENDIF()
message(STATUS "Particle properties layout (LAYOUT): " ${LAYOUT} "(-major)")

# Include third-party libraries configuration
include(TPLs)

# Make sure that the config file is in the search path
# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/Config.h.in"
  "${PROJECT_BINARY_DIR}/Base/Config.h"
)

# Remove build hostname from cmake cache
unset(HOSTNAME CACHE)

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
                    ${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/MonteCarlo
                    ${QUINOA_SOURCE_DIR}/RNG
                    ${QUINOA_SOURCE_DIR}/RNGTest
                    ${QUINOA_SOURCE_DIR}/SDE
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Paradigm
                    ${QUINOA_SOURCE_DIR}/Utilities
                    ${TPL_DIR}/include
                    ${TPL_DIR}/include/pegtl
                    ${TPL_DIR}/hdf5ser/include
                    ${PROJECT_BINARY_DIR}/Base  # for Config.h
                    ${PROJECT_BINARY_DIR}/Main  # for Charm++ modules in Main
                    ${MPI_INCLUDE_PATH}
)

# Components
add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(Geometry)
add_subdirectory(IO)
add_subdirectory(LinearAlgebra)
add_subdirectory(Mesh)
add_subdirectory(Paradigm)
add_subdirectory(MonteCarlo)
add_subdirectory(RNG)
add_subdirectory(RNGTest)
add_subdirectory(Statistics)
add_subdirectory(Main)
