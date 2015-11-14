# Build quinoa third-party Libraries

cmake_minimum_required(VERSION 2.8.5)

project(TPL C CXX Fortran)
include(ExternalProject)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
disallow_in_source_builds()

# Find MPI and get its underlying compilers
# The get_mpi_compilers cmake function uses the built-in CMake functionality to
# find the MPI wrappers and based on the wrappers detects the underyling
# compilers.
include(MPICompilers)
get_mpi_compilers()

include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(PROCESSOR_COUNT EQUAL 0)     # if ncpus could not be determined, assign 1
  set(PROCESSOR_COUNT "1")
endif()
message(STATUS "Found ${PROCESSOR_COUNT} PEs")

# Set compiler id
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)

# Set install directory based on compiler ID if was not given on command line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install/${COMP}" CACHE STRING
      "Third-party install prefix." FORCE)
  message(STATUS "CMAKE_INSTALL_PREFIX not specified, setting to '${PROJECT_SOURCE_DIR}/install/${COMP}'")
else()
  message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
endif()

# Set Charm++ build architecture: if MPI is detected, use it
if(NOT CHARM_ARCH)
  if(MPI_CXX_FOUND)
    set(CHARM_ARCHITECTURE "mpi-linux-x86_64")
    message(STATUS "CHARM_ARCH not specified and MPI found, setting to '${CHARM_ARCHITECTURE}'")
  else()
    set(CHARM_ARCHITECTURE "net-linux-x86_64")
    message(STATUS "CHARM_ARCH not specified and MPI not found, setting to '${CHARM_ARCHITECTURE}'")
  endif()

  set(CHARM_ARCH ${CHARM_ARCHITECTURE} CACHE STRING "Select Charm++ target architecture, see http://charm.cs.illinois.edu/manuals/html/charm++ and <charm>/build --help for more info.")
endif()
message(STATUS "Charm++ architecture: ${CHARM_ARCH}")

# Set Charm++ compiler (CHARM_COMPILER, the compiler used to build Charm++).
# If MPI is detected use its wrapper, if not, let the Charm++ build decide.
if(NOT CHARM_COMPILER)
  if(MPI_CXX_FOUND)
    set(CHARM_COMP "mpicxx")
    message(STATUS "CHARM_COMPILER not specified and MPI found, setting to '${CHARM_COMP}'")
  else()
    set(CHARM_COMP "")
    message(STATUS "CHARM_COMPILER not specified and MPI not found, letting Charm++ build use its default compiler")
  endif()

  set(CHARM_COMPILER ${CHARM_COMP} CACHE STRING "Select compiler used to build Charm++ with, see http://charm.cs.illinois.edu/manuals/html/charm++ and <charm>/build --help for more info.")
endif()
message(STATUS "Charm++ compiler (the compiler used to build Charm++): ${CHARM_COMPILER}")

# Compilers
# =========
# There are three sets of compilers, given by the following CMake variables:
#
# 1. The underlying compiler suite
#
#    UNDERLYING_C_COMPILER
#    UNDERLYING_CXX_COMPILER
#    UNDERLYING_Fortran_COMPILER
#
# 2. The MPI wrappers (of the underlying compilers)
#
#    MPI_C_COMPILER
#    MPI_CXX_COMPILER
#    MPI_Fortran_COMPILER
#
# 3. The Charm++ wrapper (of the MPI wrappers or around the underlying
#    compilers, depending on CHARM_COMPILER), the compiler used to build Charm++
#
#    CHARM_COMPILER

# Set Boost's default toolset
set(BOOST_TOOLSET ${COMP})

# Include function that allows to set standard C++ library
include(libstdcxx)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")  # clang-specific settings

  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

  # Don't error out on functions not returning anything even when supposed to,
  # required only by TestU01 
  set(CFLAGS "${CFLAGS} -Wno-error=return-type")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Boost as build toolset wants 'gcc' for the gnu compiler
  set(BOOST_TOOLSET "gcc")
  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libstdc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  # Link to thread-safe MPI libraries, needed for Intel MPI
  set(CFLAGS "${CFLAGS} -mt_mpi")
  set(CXXFLAGS "${CXXFLAGS} -mt_mpi")
  set(FFLAGS "${FFLAGS} -mt_mpi")
  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libstdc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")  # pgi-specific settings

  # Some configure scripts don't like pgi's preprocessor so use gnu's
  set(PGCPP "CPP=cpp")

endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'DEBUG'")
  set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Set Boost's default build type
string(TOLOWER ${CMAKE_BUILD_TYPE} BOOST_BUILD_TYPE)

# Set compiler flags for use of libc++. If nothing is given, do nothing; gnu
# and intel will use libstdc++, while clang will use libc++ by default.
if(LIBCXX_CPP)
   # Use libc++ (http://libcxx.llvm.org)
   # Use libc++abi (http://libcxxabi.llvm.org)
   set(CXXFLAGS "${CXXFLAGS} -stdlib=libc++ -lc++abi")
endif()

# Get compiler flags (CFLAGS, CXXFLAGS, FFLAGS) independent of CMAKE_BUILD_TYPE
# and echo flags that will be passed to all TPL builds
include(get_compiler_flags)
get_compiler_flags()

# Detect Boost, exclude if found, copy if not
if (NOT NO_SYSTEM_BOOST)
  find_package(Boost)
endif()
if(Boost_FOUND)
  set(boost "")
  #message(STATUS "Boost found at ${Boost_INCLUDE_DIR} - will not build ours")
else()
  set(boost "boost")
endif()

# Detect PStreams, exclude if found, build if not
if (NOT NO_SYSTEM_PSTREAMS)
  find_package(PStreams)
endif()
if(PSTREAMS_FOUND)
  #message(STATUS "PStreams found at ${PSTREAMS_INCLUDES} - will not unpack ours")
  set(pstreams "")
else()
  #message(STATUS "PStreams not found - will unpack ours")
  set(pstreams "pstreams")
endif()

# Detect Hypre, exclude if found, build if not
if (NOT NO_SYSTEM_HYPRE)
  find_package(Hypre)
endif()
if(HYPRE_FOUND)
  #message(STATUS "Hypre found at ${HYPRE_INCLUDES} (include) and at ${HYPRE_LIBRARIES} (lib) - will not build ours")
  set(hypre "")
else()
  #message(STATUS "Hypre not found - will build ours")
  set(hypre "hypre")
endif()

# Define fine-grained targets
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS
             download patch configure build install test)

#### Charm++ ###################################################################
# http://charm.cs.illinois.edu
ExternalProject_Add(
  charm
  PREFIX charm
  # Move unpacked tarball to install dir and build there
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_SOURCE_DIR}/src/charm
                    ${CMAKE_INSTALL_PREFIX}/charm
  # Build in install dir and leave stuff there - no install step
  BINARY_DIR ${CMAKE_INSTALL_PREFIX}/charm
  # Build Charm++
  BUILD_COMMAND ${CMAKE_INSTALL_PREFIX}/charm/build
                msa
                ${CHARM_ARCH}
                ${CHARM_COMPILER}
                --build-shared
                -j${PROCESSOR_COUNT}
                ${CXXFLAGS}
  DOWNLOAD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_CONFIGURE 1
  LOG_BUILD 1
)

#### PEGTL (Parsing Expression Grammar Template Library) #######################
# https://code.google.com/p/pegtl/wiki/PEGTL0
# Header only
ExternalProject_Add(
  pegtl
  PREFIX pegtl
  # Header-only, copy include dir over
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_SOURCE_DIR}/src/pegtl/include
                    ${CMAKE_INSTALL_PREFIX}/include/pegtl
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  DOWNLOAD_COMMAND ""
)

#### TUT (Template Unit Test Framework) ########################################
# http://tut-framework.sourceforge.net
# Header only
ExternalProject_Add(
  tut
  PREFIX tut
  # Header-only, copy include dir over
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_SOURCE_DIR}/src/tut/include
                    ${CMAKE_INSTALL_PREFIX}/include/tut
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  DOWNLOAD_COMMAND ""
)

#### Boost ####################################################################
# http://www.boost.org
# Header only, only if not found
if (boost)
  ExternalProject_Add(
    boost
    PREFIX boost
    # Header-only, copy include dir over
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                      ${PROJECT_SOURCE_DIR}/src/boost
                      ${CMAKE_INSTALL_PREFIX}/include/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    DOWNLOAD_COMMAND ""
  )
endif()

#### Boost.MPL extension: cartesian product ###################################
# http://www.organicvectory.com/index.php?option=com_content&view=article&id=75:boostmplcartesianproduct&catid=42:boost&Itemid=78
# Header only
ExternalProject_Add(
  cartesian_product
  PREFIX cartesian_product
  # Header-only, copy include dir over
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_SOURCE_DIR}/src/cartesian_product/boost/mpl
                    ${CMAKE_INSTALL_PREFIX}/include/boost/mpl
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  DOWNLOAD_COMMAND ""
)

#### Pstreams (POSIX Process Control in C++) ###################################
# http://pstreams.sourceforge.net/
# Header only
if (pstreams)
  ExternalProject_Add(
    pstreams
    PREFIX pstreams
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory 
                      ${CMAKE_INSTALL_PREFIX}/include/pstreams
    BUILD_COMMAND ${CMAKE_COMMAND} -E copy
                  ${PROJECT_SOURCE_DIR}/src/pstreams/pstream.h
                  ${CMAKE_INSTALL_PREFIX}/include/pstreams
    INSTALL_COMMAND ""
    DOWNLOAD_COMMAND ""
  )
endif()

#### HDF5 serial ###############################################################
# http://www.hdfgroup.org/HDF5
ExternalProject_Add(
  hdf5
  PREFIX hdf5
  DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${PROJECT_SOURCE_DIR}/src/hdf5
                   ${PROJECT_BINARY_DIR}/hdf5
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5/configure
                    CC=${MPI_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    ${PGCPP}
                    --enable-parallel
                    --enable-shared
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
) 

#### NETCDF ####################################################################
## http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
ExternalProject_Add(
  netcdf
  DEPENDS hdf5
  PREFIX netcdf
  DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${PROJECT_SOURCE_DIR}/src/netcdf
                   ${PROJECT_BINARY_DIR}/netcdf
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/netcdf/configure
                    CC=${MPI_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    ${PGCPP}
                    CPPFLAGS=-I${CMAKE_INSTALL_PREFIX}/include
                    LDFLAGS=-L${CMAKE_INSTALL_PREFIX}/lib
                    --disable-examples
                    --enable-netcdf4
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

##### Trilinos #################################################################
## http://trilinos.sandia.gov
ExternalProject_Add(
  trilinos
  DEPENDS netcdf
  PREFIX trilinos
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/trilinos
  CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
             -D TPL_ENABLE_MPI:BOOL=ON
             -D Trilinos_ENABLE_Zoltan:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodus:BOOL=ON
             -D Trilinos_ENABLE_SEACASNemesis:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodiff:BOOL=ON
             -D TPL_ENABLE_Netcdf:BOOL=ON
             -D Netcdf_LIBRARY_DIRS:PATH=${netcdf_lib}
             -D TPL_Netcdf_INCLUDE_DIRS:PATH=${netcdf_include}
             -D CMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_PREFIX}/lib
             -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_Fortran_COMPILER:PATH=${MPI_Fortran_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
             -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
             -D CMAKE_Fortran_FLAGS:STRING=${FFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  DOWNLOAD_COMMAND ""
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### Hypre #####################################################################
# http://computation.llnl.gov/project/linear_solvers/software.php
if (hypre)
  ExternalProject_Add(
    hypre
    PREFIX hypre
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/hypre/src
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_CXX_COMPILER:PATH=${UNDERLYING_CXX_COMPILER}
               -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
               -D HYPRE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    DOWNLOAD_COMMAND ""
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### RNGSSELib2 ################################################################
# http://cpc.cs.qub.ac.uk/summaries/AEIT_v2_0.html
# No need for MPI, use underlying compilers
ExternalProject_Add(
  rngsse2
  PREFIX rngsse2
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/rngsse2
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_CXX_COMPILER:PATH=${UNDERLYING_CXX_COMPILER}
             -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  DOWNLOAD_COMMAND ""
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### TestU01 ##################################################################
# http://www.iro.umontreal.ca/~simardr/testu01/tu01.html
# No need for MPI, use underlying compilers
ExternalProject_Add(
  testu01
  PREFIX testu01
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/testu01
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${UNDERLYING_C_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  DOWNLOAD_COMMAND ""
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)
