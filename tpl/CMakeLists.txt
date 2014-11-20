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

# Set Charm++/AMPI build architecture: if MPI is detected, use it
if(NOT CHARM_ARCH)
  if(MPI_CXX_FOUND)
    set(CHARM_ARCHITECTURE "mpi-linux-x86_64")
    message(STATUS "CHARM_ARCH not specified and MPI found, setting to '${CHARM_ARCHITECTURE}'")
  else()
    set(CHARM_ARCHITECTURE "net-linux-x86_64")
    message(STATUS "CHARM_ARCH not specified and MPI not found, setting to '${CHARM_ARCHITECTURE}'")
  endif()

  set(CHARM_ARCH ${CHARM_ARCHITECTURE} CACHE STRING "Select Charm++/AMPI target architecture, see http://charm.cs.illinois.edu/manuals/html/charm++ and <charm>/build --help for more info.")
endif()

# Set Charm++/AMPI compiler (the compiler used to build Charm++/AMPI): if MPI
# is detected use it
if(NOT CHARM_COMPILER)
  if(MPI_CXX_FOUND)
    set(CHARM_COMP "mpicxx")
    message(STATUS "CHARM_COMPILER not specified and MPI found, setting to '${CHARM_COMP}'")
  else()
    set(CHARM_COMP "")
    message(STATUS "CHARM_COMPILER not specified and MPI not found, letting Charm++ build use its default compiler")
  endif()

  set(CHARM_COMPILER ${CHARM_COMP} CACHE STRING "Select compiler used to build Charm++/AMPI with, see http://charm.cs.illinois.edu/manuals/html/charm++ and <charm>/build --help for more info.")
endif()

# Set Charm++/AMPI compiler wrappers
set(CHARM_AMPI_C_COMPILER ${CMAKE_INSTALL_PREFIX}/charm/bin/ampicc)
set(CHARM_AMPI_CXX_COMPILER ${CMAKE_INSTALL_PREFIX}/charm/bin/ampicxx)
set(CHARM_AMPI_Fortran_COMPILER ${CMAKE_INSTALL_PREFIX}/charm/bin/ampif90)

# Echo Charm++ underlying compilers
MESSAGE(STATUS "Charm++/AMPI C compiler: " ${CHARM_AMPI_C_COMPILER})
MESSAGE(STATUS "Charm++/AMPI C++ compiler: " ${CHARM_AMPI_CXX_COMPILER})
MESSAGE(STATUS "Charm++/AMPI Fortran compiler: " ${CHARM_AMPI_Fortran_COMPILER})

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
# 3. The Charm++/AMPI wrappers (of the MPI wrappers or around the underlying
#    compilers, depending on CHARM_COMPILER), the compiler used to build Charm++
#    and AMPI
#
#    CHARM_AMPI_C_COMPILER
#    CHARM_AMPI_CXX_COMPILER
#    CHARM_AMPI_Fortran_COMPILER

# Set Boost's default toolset
set(BOOST_TOOLSET ${COMP})

# Include function that allows to set standard C++ library
include(libstdcxx)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")  # clang-specific settings

  # Set libc++ as the default standard C++ library and allow selection
  setlibstdcxx("libc++")
  # Don't error out on functions not returning anything even when supposed to,
  # required only by TestU01 
  set(CFLAGS "${CFLAGS} -Wno-error=return-type")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Boost as build toolset wants 'gcc' for the gnu compiler
  set(BOOST_TOOLSET "gcc")
  # Set libstdc++ as the default standard C++ library and allow selection
  setlibstdcxx("libstdc++")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  # Link to thread-safe MPI libraries, needed for Intel MPI
  set(CFLAGS "${CFLAGS} -mt_mpi")
  set(CXXFLAGS "${CXXFLAGS} -mt_mpi")
  set(FFLAGS "${FFLAGS} -mt_mpi")
  # Set libstdc++ as the default standard C++ library and allow selection
  setlibstdcxx("libstdc++")

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

# Get compiler flags (CFLAGS, CXXFLAGS, FFLAGS) independent of CMAKE_BUILD_TYPE
include(get_compiler_flags)
get_compiler_flags()

# Set compiler flags for use of libc++ (If nothing is given, do nothing; gnu,
# intel, and clang will use libstdc++ by default).
if(LIBCXX_CPP)
  set(CXXFLAGS_LIBCPP "-stdlib=libc++")
  set(CXXFLAGS_LIBCPPABI "-lc++abi")
  set(CXXFLAGS "${CXXFLAGS} ${CXXFLAGS_LIBCPP} ${CXXFLAGS_LIBCPPABI}")
endif()

# Detect Boost, exclude if found, copy if not
find_package(Boost)
if(Boost_FOUND)
  set(boost "")
  message(STATUS "Boost found at ${Boost_INCLUDE_DIR} - will not build ours")
else()
  set(boost "boost")
endif()

# Detect HDF5, exclude if found, build if not
find_package(HDF5)
if(HDF5_FOUND)
  message(STATUS "HDF5 found at ${HDF5_INCLUDE_DIRS} (include) and at ${HDF5_LIBRARY_DIRS} (libs) - will not build ours")
  set(hdf5ser "")
  set(hdf5ser_include "${HDF5_INCLUDE_DIRS}")
  set(hdf5ser_lib "${HDF5_LIBRARY_DIRS}")
else()
  set(hdf5ser "hdf5ser")
  set(hdf5ser_include "${CMAKE_INSTALL_PREFIX}/hdf5ser/include")
  set(hdf5ser_lib "${CMAKE_INSTALL_PREFIX}/hdf5ser/lib")
endif()

# Define fine-grained targets
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS
             download patch configure build install test)

#### Charm++, AMPI #############################################################
# http://charm.cs.illinois.edu
ExternalProject_Add(
  charm
  PREFIX charm
  URL ${PROJECT_SOURCE_DIR}/charm-6.6.0.tar.gz
  URL_MD5 31e95901b3f7324d52107e6ad000fcc8
  # Move unpacked tarball to install dir and build there
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_BINARY_DIR}/charm/src/charm
                    ${CMAKE_INSTALL_PREFIX}/charm
  # Build in install dir and leave stuff there - no install step
  BINARY_DIR ${CMAKE_INSTALL_PREFIX}/charm
  # Build AMPI which triggers building Charm++
  BUILD_COMMAND ${CMAKE_INSTALL_PREFIX}/charm/build
                AMPI
                ${CHARM_ARCH}
                ${CHARM_COMPILER}
                --build-shared
                -j8
                ${CXXFLAGS}
                ${CXXFLAGS_LIBCPP} ${CXXFLAGS_LIBCPPABI}
  INSTALL_COMMAND ""
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
)

#### PEGTL (Parsing Expression Grammar Template Library) #######################
# http://code.google.com/p/pegtl
# Header only
ExternalProject_Add(
  pegtl
  PREFIX pegtl
  URL ${PROJECT_SOURCE_DIR}/pegtl-0.32.tar.gz
  URL_MD5 b2989e2b38f0a80f197f59ed00c50a81
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD 1
)
# Header-only, copy include dir over
ExternalProject_Get_Property(pegtl SOURCE_DIR)
ExternalProject_Add_Step(
  pegtl copy-source
  DEPENDEES download
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${SOURCE_DIR}/include
          ${CMAKE_INSTALL_PREFIX}/include/pegtl
)

#### TUT (Template Unit Test Framework) ########################################
# http:://tut-framework.sourceforge.net
# Header only
ExternalProject_Add(
  tut
  PREFIX tut
  URL ${PROJECT_SOURCE_DIR}/tut-framework-2013-12-18.tar.gz
  URL_MD5 49ac872ee8d6066fd5ee499b38adbb11
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD 1
)
# Header-only, copy include dir over
ExternalProject_Get_Property(tut SOURCE_DIR)
ExternalProject_Add_Step(
  tut copy-source
  DEPENDEES download
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${SOURCE_DIR}/include
          ${CMAKE_INSTALL_PREFIX}/include/tut
)

#### Boost ####################################################################
# http://www.boost.org
# Header only, only if not found
if (boost)
  ExternalProject_Add(
    boost
    PREFIX boost
    URL ${PROJECT_SOURCE_DIR}/boost_1_54_0.tar.bz2
    URL_MD5 15cb8c0803064faef0c4ddf5bc5ca279
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 1
  )
  # Header-only, copy include dir over
  ExternalProject_Get_Property(boost SOURCE_DIR)
  ExternalProject_Add_Step(
    boost copy-source
    DEPENDEES download
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${SOURCE_DIR}/boost
            ${CMAKE_INSTALL_PREFIX}/include/boost
  )
endif()

#### Boost.MPL extension: cartesian product ###################################
# http://www.organicvectory.com/index.php?option=com_content&view=article&id=75:boostmplcartesianproduct&catid=42:boost&Itemid=78
# Header only
ExternalProject_Add(
  cartesian_product
  PREFIX cartesian_product
  URL ${PROJECT_SOURCE_DIR}/boost.mpl.cartesian_product.zip
  URL_MD5 4bcfbdd61af85548e684198c0ee52b62
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD 1
)
# Header-only, copy include dir over
ExternalProject_Get_Property(cartesian_product SOURCE_DIR)
ExternalProject_Add_Step(
  cartesian_product copy-source
  DEPENDEES download
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${SOURCE_DIR}/boost/mpl
          ${CMAKE_INSTALL_PREFIX}/include/boost/mpl
)

#### HDF5 serial ###############################################################
# http://www.hdfgroup.org/HDF5
if (hdf5ser)
  ExternalProject_Add(
    hdf5ser
    PREFIX hdf5ser
    URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.13.tar.gz
    URL_MD5 c03426e9e77d7766944654280b467289
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5ser/src/hdf5ser/configure
                      CC=${UNDERLYING_C_COMPILER}
                      CFLAGS=${CFLAGS}
                      ${PGCPP}
                      --prefix=${CMAKE_INSTALL_PREFIX}/hdf5ser
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
  )
endif()

#### HDF5 parallel #############################################################
# http://www.hdfgroup.org/HDF5
#ExternalProject_Add(
#  hdf5par
#  PREFIX hdf5par
#  URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.13.tar.gz
#  URL_MD5 c03426e9e77d7766944654280b467289
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5par/src/hdf5par/configure
#                    CC=${MPI_C_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    --enable-parallel
#                    --prefix=${CMAKE_INSTALL_PREFIX}/hdf5par
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

##### H5Part ###################################################################
## http://vis.lbl.gov/Research/H5Part
#ExternalProject_Add(
#  h5part
#  DEPENDS hdf5ser
#  PREFIX h5part
#  URL ${PROJECT_SOURCE_DIR}/H5Part-1.6.6.tar.gz
#  URL_MD5 327c63d198e38a12565b74cffdf1f9d7
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/h5part/src/h5part/configure
#                    CC=${MPI_C_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    --enable-parallel
#                    --with-hdf5=${CMAKE_INSTALL_PREFIX}/hdf5ser
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### zlib ######################################################################
# http://www.zlib.net
#ExternalProject_Add(
#  zlib
#  PREFIX zlib
#  URL ${PROJECT_SOURCE_DIR}/zlib-1.2.8.tar.gz
#  URL_MD5 44d667c142d7cda120332623eab69f40
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/zlib/src/zlib/configure
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_COMMAND make CC=${UNDERLYING_C_COMPILER} CFLAGS=${CFLAGS} install
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### NETCDF ####################################################################
## http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
ExternalProject_Add(
  netcdf
  DEPENDS ${hdf5ser}
  PREFIX netcdf
  URL ${PROJECT_SOURCE_DIR}/netcdf-4.3.2.tar.gz
  URL_MD5 2fd2365e1fe9685368cd6ab0ada532a0
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/netcdf/src/netcdf/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    ${PGCPP}
                    CPPFLAGS=-I${hdf5ser_include}
                    LDFLAGS=-L${hdf5ser_lib}
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
  URL ${PROJECT_SOURCE_DIR}/trilinos-11.6.1-Source.tar.gz
  URL_MD5 e933ca40a585c7ba9edfc0f45186aaf7
  CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
             -D Trilinos_ENABLE_Fortran:BOOL=OFF
             -D Trilinos_ENABLE_Zoltan:BOOL=ON
             -D TPL_ENABLE_MPI:BOOL=ON
             -D Trilinos_ENABLE_SEACASExodus:BOOL=ON
             -D Trilinos_ENABLE_SEACASNemesis:BOOL=ON
             -D TPL_ENABLE_Netcdf:BOOL=ON
             -D Netcdf_LIBRARY_DIRS:PATH=${CMAKE_INSTALL_PREFIX}/lib
             -D TPL_Netcdf_INCLUDE_DIRS:PATH=${CMAKE_INSTALL_PREFIX}/include
             -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
             -D CMAKE_CXX_FLAGS:STRING=${CXXFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

##### Silo #####################################################################
## https://wci.llnl.gov/codes/silo
#ExternalProject_Add(
#  silo
#  DEPENDS hdf5ser
#  PREFIX silo
#  URL ${PROJECT_SOURCE_DIR}/silo-4.9.1-bsd.tar.gz
#  URL_MD5 a66223e1c0e569eb4887df2aae0ab70f
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/silo/src/silo/configure
#                    CC=${MPI_C_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    ${PGCPP}
#                    --enable-optimization
#                    --enable-shared
#                    --disable-silex
#                    --with-hdf5=${CMAKE_INSTALL_PREFIX}/hdf5ser/include,${CMAKE_INSTALL_PREFIX}/hdf5ser/lib
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### minixml ###################################################################
# http://www.minixml.org/software.php
#ExternalProject_Add(
#  minixml
#  PREFIX minixml
#  URL ${PROJECT_SOURCE_DIR}/mxml-2.8.tar.gz
#  URL_MD5 d85ee6d30de053581242c4a86e79a5d2
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/minixml/src/minixml/configure
#                    CC=${UNDERLYING_C_COMPILER}
#                    CXX=${UNDERLYING_CXX_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    CXXFLAGS=${CXXFLAGS}
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### bzip2 #####################################################################
# http://www.bzip.org
#ExternalProject_Add(
#  bzip2
#  PREFIX bzip2
#  URL ${PROJECT_SOURCE_DIR}/bzip2-1.0.6.tar.gz
#  URL_MD5 00b516f4704d4a7cb50a1d97e6e8e15b
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND make
#                CC=${UNDERLYING_C_COMPILER}
#                CFLAGS=${CFLAGS}
#                PREFIX=${CMAKE_INSTALL_PREFIX}
#                install
#  INSTALL_COMMAND ""
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### szip ######################################################################
# http://www.hdfgroup.org/ftp/lib-external/szip
#ExternalProject_Add(
#  szip
#  PREFIX szip
#  URL ${PROJECT_SOURCE_DIR}/szip-2.1.tar.gz
#  URL_MD5 902f831bcefb69c6b635374424acbead
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/szip/src/szip/configure
#                    CC=${UNDERLYING_C_COMPILER}
#                    CXX=${UNDERLYING_CXX_COMPILER}
#                    FC=${UNDERLYING_Fortran_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    CXXFLAGS=${CXXFLAGS}
#                    FFLAGS=${FFLAGS}
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### ISOBAR ####################################################################
# http://dx.doi.org/10.1145/2287076.2287086
#ExternalProject_Add(
#  isobar
#  PREFIX isobar
#  URL ${PROJECT_SOURCE_DIR}/isobar.0.3.0.tgz
#  URL_MD5 86d123b7028ed238c0e6e49c22e7272c
#  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
#            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include
#  BUILD_COMMAND make install
#  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
#                  ${PROJECT_BINARY_DIR}/isobar/src/isobar/lib/libisobar.a
#                  ${CMAKE_INSTALL_PREFIX}/lib
#          COMMAND ${CMAKE_COMMAND} -E copy
#                  ${PROJECT_BINARY_DIR}/isobar/src/isobar/include/isobar.h
#                  ${CMAKE_INSTALL_PREFIX}/include
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### APLOD #####################################################################
# http://conferences.computer.org/sc/2012/papers/1000a088.pdf
#ExternalProject_Add(
#  aplod
#  PREFIX aplod
#  URL ${PROJECT_SOURCE_DIR}/aplod.0.1.tar.gz
#  URL_MD5 3e9749efb42a3d4d3165e884456812bc
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND make CXX=${UNDERLYING_CXX_COMPILER} CXXLAGS=${CXXFLAGS}
#  INSTALL_COMMAND make install prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### ADIOS #####################################################################
# https://www.olcf.ornl.gov/center-projects/adios
#ExternalProject_Add(
#  adios
#  DEPENDS minixml zlib bzip2 szip isobar aplod hdf5ser hdf5par netcdf
#  PREFIX adios
#  URL ${PROJECT_SOURCE_DIR}/adios-1.6.0.tar.gz
#  URL_MD5 70073584c620b92d7842250ce778847f
#  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/adios/src/adios/configure
#                    CC=${MPI_C_COMPILER}
#                    CXX=${MPI_CXX_COMPILER}
#                    CFLAGS=${CFLAGS}
#                    CXXLAGS=${CXXFLAGS}
#                    --disable-fortran
#                    --enable-research-transports
#                    --enable-skel-timing
#                    --with-mxml=${CMAKE_INSTALL_PREFIX}
#                    --with-phdf5=${CMAKE_INSTALL_PREFIX}/hdf5par
#                    --with-hdf5=${CMAKE_INSTALL_PREFIX}/hdf5ser
#                    --with-netcdf=${CMAKE_INSTALL_PREFIX}
#                    --with-zlib=${CMAKE_INSTALL_PREFIX}
#                    --with-bzip2=${CMAKE_INSTALL_PREFIX}
#                    --with-szip=${CMAKE_INSTALL_PREFIX}
#                    --with-isobar=${CMAKE_INSTALL_PREFIX}
#                    --with-aplod=${CMAKE_INSTALL_PREFIX}
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)

#### RNGSSELib2 ################################################################
# http://dx.doi.org/10.1016/j.cpc.2013.04.007
# No need for MPI, use underlying compilers
ExternalProject_Add(
  rngsse2
  PREFIX rngsse2
  URL ${PROJECT_SOURCE_DIR}/rngsselib_2.tar.gz
  URL_MD5 dc0db58196dd291a362028545b55e2a9
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_CXX_COMPILER:PATH=${UNDERLYING_CXX_COMPILER}
             -D CMAKE_CXX_FLAGS:STRING=${CXFLAGS}
             -D CMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  LOG_DOWNLOAD 1
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
  URL ${PROJECT_SOURCE_DIR}/TestU01-1.2.3.tar.gz
  URL_MD5 dfe63dfce27d1a6114a644f6eb9a23a7
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/testu01/src/testu01/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CXX=${UNDERLYING_CXX_COMPILER}
                    F77=${UNDERLYING_Fortran_COMPILER}
                    ${PGCPP}
                    CFLAGS=${CFLAGS}
                    CPPFLAGS=${CFLAGS}
                    CXXFLAGS=${CXXFLAGS}
                    FFLAGS=${FFLAGS}
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### TestU01-hack ##############################################################
# http://www.iro.umontreal.ca/~simardr/testu01/tu01.html
# No need for MPI, use underlying compilers
# This TPL-build can be used to hack TestU01. ExternalProject_Add will look for
# the source in SOURCE_DIR and build there. Building with different compilers
# requires a 'make clean' before this cmake command will be successful. No
# PREFIX, no URL, no URL_MD5, instead we specify SOURCE_DIR, containing the
# fork. Changes are tracked in git repo in SOURCE_DIR. Once the hacking is done,
# create a new tarball and update the above ExternalProject_Add, as that is more
# automatic for new builds (using different compilers) from scratch.
#ExternalProject_Add(
#  testu01
#  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/TestU01-1.2.3
#  CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/src/TestU01-1.2.3/configure
#                    CC=${UNDERLYING_C_COMPILER}
#                    CXX=${UNDERLYING_CXX_COMPILER}
#                    F77=${UNDERLYING_Fortran_COMPILER}
#                    CFLAGS=-Wno-error=return-type
#                    CPPFLAGS=${CFLAGS}
#                    CXXFLAGS=${CXXFLAGS}
#                    FFLAGS=${FFLAGS}
#                    --prefix=${CMAKE_INSTALL_PREFIX}
#  BUILD_IN_SOURCE 1
#  LOG_DOWNLOAD 1
#  LOG_CONFIGURE 1
#  LOG_BUILD 1
#  LOG_INSTALL 1
#)
