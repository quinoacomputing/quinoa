################################################################################
#
# \file      src/CMakeLists.txt
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Build quinoa third-party Libraries
# \date      Fri 06 May 2016 06:46:49 AM MDT
#
################################################################################

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(TPL C CXX Fortran)
include(ExternalProject)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
disallow_in_source_builds()

# Set CMAKE_BUILD_TYPE
include(BuildType)

# Set BUILD_SHARED_LIBS
include(BuildShared)

# Conditionally enforce static linking for exeuctables added after this point.
# This must precede all add_executable() calls adding executables to link
# statically.
if (NOT BUILD_SHARED_LIBS)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

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

# Include support for multiarch path names
include(GNUInstallDirs)

# Set compiler id
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)

# Detect operating system type and version
include(DetectOS)
detect_os()
message(STATUS "Detected OS: ${HOST_OS}")

# Query target architecture
include(TargetArch)
target_architecture(ARCH)
message(STATUS "Detected target architecture: ${ARCH}")

# Set install directory based on compiler ID if was not given on command line
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install/${COMP}-${ARCH}"
      CACHE STRING "Third-party install prefix." FORCE)
  message(STATUS "CMAKE_INSTALL_PREFIX not specified, setting to '${PROJECT_SOURCE_DIR}/install/${COMP}-${ARCH}'")
else()
  message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
endif()

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

# Convert system name to lower case
string(TOLOWER ${CMAKE_SYSTEM_NAME} SYSTEM_NAME)

# Match Charm++ build architecture for system architecture detected
if(ARCH MATCHES "x86_64")

  set(CHARM_ARCH "-${ARCH}")

elseif(ARCH MATCHES "i386")

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CHARM_ARCH "")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CHARM_ARCH "-x86")
  endif()

else()

  message(FATAL "Selecting a Charm++ architecture for ${ARCH} is not supported")

endif()

# Set Charm++ architecture
set(CHARM_ARCHITECTURE "mpi-${SYSTEM_NAME}${CHARM_ARCH}")

# Set Charm++ compiler
if(CMAKE_CXX_COMPILER MATCHES "mpi")    # using mpi wrappers

  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CHARM_COMPILER "mpicxx")
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CHARM_COMPILER "")
  endif()

else()  # if not using mpi wrappers

  set(CHARM_COMPILER "")

endif()

message(STATUS "Charm++ architecture: '${CHARM_ARCHITECTURE}'")
message(STATUS "Charm++ compiler: '${CHARM_COMPILER}'")

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")  # clang-specific settings

  # Include function that allows to set standard C++ library
  include(libstdcxx)

  # Offer switch between libc++ and libstdc++, using libc++ as default
  if (libc++_FOUND)
    set_libstdcpp_vs_libcpp(CXXFLAGS "libc++")
  endif()

  # Don't error out on functions not returning anything even when supposed to,
  # required only by TestU01 
  set(CFLAGS "${CFLAGS} -Wno-error=return-type")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Gnu will use libstdc++ by default. We do not offer libc++.

  # Detect gcc version
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
  list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
  list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
  message(STATUS "GCC version: ${GCC_MAJOR}.${GCC_MINOR}")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  # Include function that allows to set standard C++ library
  include(libstdcxx)

  # Offer switch between libc++ and libstdc++, using libc++ as default
  if (libc++_FOUND)
    set_libstdcpp_vs_libcpp(CXXFLAGS "libc++")
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")  # pgi-specific settings

  # Some configure scripts don't like pgi's preprocessor so use gnu's
  set(PGCPP "CPP=cpp")

endif()

# Get compiler flags (CFLAGS, CXXFLAGS, FFLAGS) independent of CMAKE_BUILD_TYPE
# and echo flags that will be passed to all TPL builds
include(get_compiler_flags)
get_compiler_flags()

# Mac OS X specifics
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX_RPATH "-DCMAKE_MACOSX_RPATH:BOOL=ON")
endif()

# Attempt to find some TPLs on the system (some others we always build)

# Detect HDF5, exclude if found, copy if not
set(HDF5_PREFER_PARALLEL true)
if(NOT NO_SYSTEM_HDF5)
  if(NOT BUILD_SHARED_LIBS)
    set(HDF5_USE_STATIC_LIBRARIES true)
  endif()
  find_package(HDF5 COMPONENTS HL)
endif()
if(HDF5_FOUND AND HDF5_IS_PARALLEL)
  set(hdf5 "")
  set(hdf5_include_dirs_flags "${HDF5_INCLUDE_DIRS}")
  string(REGEX REPLACE "([^;]+)" "-I\\1" hdf5_include_dirs_flags
         "${hdf5_include_dirs_flags}")
  set(hdf5_library_dirs_flags "${HDF5_LIBRARY_DIRS}")
  string(REGEX REPLACE "([^;]+)" "-L\\1" hdf5_library_dirs_flags
        "${hdf5_library_dirs_flags}")
  set(hdf5_library_dirs "${HDF5_LIBRARY_DIRS}")
else()
  message(STATUS "HDF5 either not found or not parallel: will build")
  set(hdf5 "hdf5")
  set(hdf5_include_dirs_flags "-I${CMAKE_INSTALL_PREFIX}/include")
  set(hdf5_library_dirs_flags "-L${CMAKE_INSTALL_PREFIX}/lib")
  set(hdf5_library_dirs "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Detect NetCDF, exclude if found, build if not
if(NOT NO_SYSTEM_NETCDF)
  find_package(NetCDF)
endif()
if(NETCDF_FOUND)
  set(netcdf "")
  set(netcdf_include_dirs "${NETCDF_INCLUDES}")
  set(netcdf_library_dirs "${NETCDF_LIBRARY_DIRS}")
else()
  set(netcdf "netcdf")
  set(netcdf_include_dirs "${CMAKE_INSTALL_PREFIX}/include")
  set(netcdf_library_dirs "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Detect Boost, exclude if found, copy if not
if(NOT NO_SYSTEM_BOOST)
  if(NOT BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS true)
  endif()
  find_package(Boost)
endif()
if(Boost_FOUND)
  set(boost "")
else()
  set(boost "boost")
endif()

# Detect PStreams, exclude if found, build if not
if(NOT NO_SYSTEM_PSTREAMS)
  find_package(PStreams)
endif()
if(PSTREAMS_FOUND)
  set(pstreams "")
else()
  set(pstreams "pstreams")
endif()

# Detect Hypre, exclude if found, build if not
if(NOT NO_SYSTEM_HYPRE)
  find_package(Hypre)
endif()
if(HYPRE_FOUND)
  set(hypre "")
else()
  set(hypre "hypre")
endif()

# Detect pugixml, exclude if found, build if not
if(NOT NO_SYSTEM_PUGIXML)
  find_package(pugixml)
endif()
if(PUGIXML_FOUND)
  set(pugixml "")
else()
  set(pugixml "pugixml")
endif()

# Detect blas and lapack, exclude if found, build if not
set(lapack "lapack")
set(blas "blas")
find_package(MKL)     # Prefer MKL if found
if(MKL_FOUND)
  set(blas "")
  set(lapack "")
else()
  if(NOT NO_SYSTEM_BLAS)
    find_package(BLAS)
    if(BLAS_FOUND)
      set(blas "")
    endif()
  endif()
  if(NOT NO_SYSTEM_LAPACK)
    find_package(LAPACKE)
    if(LAPACKE_FOUND)
      set(lapack "")
    endif()
  endif()
endif()

# Detect AEC (only for static link)
if(NOT BUILD_SHARED_LIBS)
  set(aec "aec")
  if(NOT NO_SYSTEM_AEC)
    find_package(AEC)
    if(AEC_FOUND)
      set(aec "")
   endif()
  endif()
endif()

# Output TPLs to be built
message(STATUS "Will build: charm trilinos pegtl tut cartesian_product rngsse2 testu01 ${hdf5} ${netcdf} ${aec} ${blas} ${boost} ${pstreams} ${hypre} ${pugixml} ${lapack}")

# Define fine-grained targets
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS
             download patch configure build install test)

#### Charm++ ###################################################################
# http://charm.cs.illinois.edu
if(BUILD_SHARED_LIBS)
  set(CHARM_SHARED "--build-shared")
endif()
# Enable "--with-production" for Charm++ build for Release builds
if(CMAKE_BUILD_TYPE MATCHES "REL")      # Release | RelWithDebInfo | MinSizeRel
  set(CHARM_RELEASE_BUILD "--with-production")
endif()
if (HOST_OS MATCHES "alpine")
  # To avoid requiring rpc.h on Alpine Linux using musl-libc, see also
  # https://lists.cs.illinois.edu/lists/arc/charm/2016-04/msg00002.html
  set(CHARM_ALPINE "--disable-charmdebug")
  # Use pthreads and avoid requiring glibc's thread context API, unimplemented
  # on purpose on Alpine Linux and musl-libc, see also the thread at
  # http://www.openwall.com/lists/musl/2016/02/04/3
  set(CHARM_PTHREADS "pthreads")
endif()
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
                charm++
                ${CHARM_ARCHITECTURE}
                ${CHARM_COMPILER}
                ${CHARM_PTHREADS}
                ${CHARM_SHARED}
                ${CHARM_RELEASE_BUILD}
                ${CHARM_ALPINE}
                -j${PROCESSOR_COUNT}
                ${CXXFLAGS}
  DOWNLOAD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_CONFIGURE 1
  LOG_BUILD 1
)

#### HDF5 ######################################################################
# http://www.hdfgroup.org/HDF5
if (hdf5)
  if(BUILD_SHARED_LIBS)
    set(HDF5_SHARED "--enable-shared")
  else()
    set(HDF5_SHARED "--disable-shared;--enable-static;--with-pic")
  endif()
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
                      ${HDF5_SHARED}
                      --prefix=${CMAKE_INSTALL_PREFIX}
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### NETCDF ####################################################################
## http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
if (netcdf)
  if(BUILD_SHARED_LIBS)
    set(NETCDF_SHARED "--enable-shared")
  else()
    set(NETCDF_SHARED "--enable-static;--disable-shared;--with-pic")
  endif()
  ExternalProject_Add(
    netcdf
    DEPENDS ${hdf5}
    PREFIX netcdf
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory
                     ${PROJECT_SOURCE_DIR}/src/netcdf
                     ${PROJECT_BINARY_DIR}/netcdf
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/netcdf/configure
                      CC=${MPI_C_COMPILER}
                      CFLAGS=${CFLAGS}
                      ${PGCPP}
                      CPPFLAGS=${hdf5_include_dirs_flags}
                      LDFLAGS=${hdf5_library_dirs_flags}
                      --enable-netcdf4
                      --disable-examples
                      --disable-fsync
                      --disable-dap
                      ${NETCDF_SHARED}
                      --prefix=${CMAKE_INSTALL_PREFIX}
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### AEC #######################################################################
## https://gitlab.dkrz.de/k202009/libaec
if (aec)
  if(BUILD_SHARED_LIBS)
    set(AEC_SHARED "--enable-shared")
  else()
    set(AEC_SHARED "--enable-static;--disable-shared")
  endif()
  ExternalProject_Add(
    aec
    PREFIX aec
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/aec
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_Fortran_COMPILER:PATH=${UNDERLYING_Fortran_COMPILER}
               -D CMAKE_Fortran_FLAGS:STRING=${FFLAGS}
               -D BUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    DOWNLOAD_COMMAND ""
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### blas ######################################################################
# http://www.netlib.org/blas
if (blas)
  ExternalProject_Add(
    blas
    PREFIX blas
    SOURCE_DIR ${PROJECT_BINARY_DIR}/blas
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory
                     ${PROJECT_SOURCE_DIR}/src/blas
                     ${PROJECT_BINARY_DIR}/blas
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory 
                      ${CMAKE_INSTALL_PREFIX}/lib
    BUILD_COMMAND make -j${PROCESSOR_COUNT} PLAT=''
                  FORTRAN=${UNDERLYING_Fortran_COMPILER}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
                    ${PROJECT_BINARY_DIR}/blas/blas.a
                    ${CMAKE_INSTALL_PREFIX}/lib/libblas.a
    BUILD_IN_SOURCE 1
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### lapack ####################################################################
# http://www.netlib.org/lapack
if (lapack)
  ExternalProject_Add(
    lapack
    PREFIX lapack
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/lapack
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_Fortran_COMPILER:PATH=${UNDERLYING_Fortran_COMPILER}
               -D CMAKE_Fortran_FLAGS:STRING=${FFLAGS}
               -D BUILD_TESTING:BOOL=OFF
               -D BUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
               -D LAPACKE:BOOL=ON
               -D USE_OPTIMIZED_BLAS:BOOL=ON
               -D USE_OPTIMIZED_LAPACK:BOOL=ON
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
               ${MACOSX_RPATH}
    DOWNLOAD_COMMAND ""
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

##### Trilinos #################################################################
## http://trilinos.sandia.gov
if(BUILD_SHARED_LIBS)
  set(TRILINOS_SHARED "-DBUILD_SHARED_LIBS=ON")
else()
  set(TRILINOS_SHARED "-DBUILD_SHARED_LIBS=OFF;-DTPL_FIND_SHARED_LIBS:BOOL=OFF;-DTrilinos_LINK_SEARCH_START_STATIC:BOOL=ON")
  set(tr_hdf5 "-DTrilinos_EXTRA_LINK_FLAGS=-L${hdf5_library_dirs} -lhdf5_hl -lhdf5 -lz -lm -ldl")
endif()
if (HOST_OS MATCHES "alpine")
  set(TRILINOS_ALPINE "-DBLAS_LIBRARY_DIRS=${CMAKE_INSTALL_PREFIX}/lib;-DLAPACK_LIBRARY_DIRS=${CMAKE_INSTALL_PREFIX}/lib")
endif()
ExternalProject_Add(
  trilinos
  DEPENDS ${hdf5} ${netcdf} ${blas} ${lapack}
  PREFIX trilinos
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/trilinos
  CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
             -D TPL_ENABLE_MPI:BOOL=ON
             -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF
             -D Trilinos_ENABLE_Zoltan2:BOOL=ON
             -D Zoltan2_ENABLE_Experimental:BOOL=ON
             -D Trilinos_ENABLE_SECONDARY_STABLE_CODE:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodus:BOOL=ON
             -D Trilinos_ENABLE_SEACASNemesis:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodiff:BOOL=ON
             -D Trilinos_ENABLE_Fortran=OFF
             -D Trilinos_ENABLE_DEBUG=OFF
             ${tr_hdf5}
             ${TRILINOS_SHARED}
             ${TRILINOS_ALPINE}
             ${MACOSX_RPATH}
             -D TPL_ENABLE_Netcdf:BOOL=ON
             -D Netcdf_LIBRARY_DIRS:PATH=${netcdf_library_dirs}
             -D TPL_Netcdf_INCLUDE_DIRS:PATH=${netcdf_include_dirs}
             -D CMAKE_INSTALL_RPATH:PATH=${CMAKE_INSTALL_PREFIX}/lib
             -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
             -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  DOWNLOAD_COMMAND ""
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
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

#### RNGSSELib2 ################################################################
# http://cpc.cs.qub.ac.uk/summaries/AEIT_v2_0.html
# No need for MPI, use underlying compilers
ExternalProject_Add(
  rngsse2
  PREFIX rngsse2
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/rngsse2
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
             ${MACOSX_RPATH}
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
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")  # force static libs on MAC OS X
  set(testu01_shared_off_mac "-DBUILD_SHARED_LIBS:BOOL=OFF")
endif()
ExternalProject_Add(
  testu01
  PREFIX testu01
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/testu01
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
             ${testu01_shared_off_mac}
             ${MACOSX_RPATH}
             -D CMAKE_C_COMPILER:PATH=${UNDERLYING_C_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  DOWNLOAD_COMMAND ""
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
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

#### Hypre #####################################################################
# http://computation.llnl.gov/project/linear_solvers/software.php
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")  # force static libs on MAC OS X
  set(hypre_shared_off_mac "-DHYPRE_SHARED:BOOL=OFF")
endif()
if (hypre)
  ExternalProject_Add(
    hypre
    PREFIX hypre
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/hypre/src
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
               -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
               ${MACOSX_RPATH}
               -D HYPRE_SHARED:BOOL=${BUILD_SHARED_LIBS}
               ${hypre_shared_off_mac}
               -D HYPRE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    DOWNLOAD_COMMAND ""
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### pugixml ###################################################################
# http://pugixml.org
if (pugixml)
  ExternalProject_Add(
    pugixml
    PREFIX pugixml
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/pugixml/scripts
    CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_CXX_COMPILER:PATH=${UNDERLYING_CXX_COMPILER}
               -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
               -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
               -D BUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
               ${MACOSX_RPATH}
    DOWNLOAD_COMMAND ""
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()
