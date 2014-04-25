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
include(GetMPICompilers)
get_mpi_compilers()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'DEBUG'")
  set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

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

# Set Boost's default toolset
set(BOOST_TOOLSET ${COMP})

# Set Boost's default build type
string(TOLOWER ${CMAKE_BUILD_TYPE} BOOST_BUILD_TYPE)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")  # clang-specific settings

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings
  # Boost as build toolset wants 'gcc' for the gnu compiler
  set(BOOST_TOOLSET "gcc")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings
  # Link to thread-safe MPI libraries
  set(CFLAGS "${CFLAGS} -mt_mpi")
  set(CXXFLAGS "${CXXFLAGS} -mt_mpi")
  set(FFLAGS "${FFLAGS} -mt_mpi")

endif()

# Convert CMAKE_DEBUG_TYPE to upper case so even if the user gives it via the
# command line, the tests below work as intended
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# Forward cmake compiler flags to automake builds
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CFLAGS "${CFLAGS} ${CMAKE_C_FLAGS}")
  set(CXXFLAGS "${CXXFLAGS} ${CMAKE_CXX_FLAGS}")
  set(FFLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CFLAGS "${CFLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(CXXFLAGS "${CXXFLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(FFLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  set(CFLAGS "${CFLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CXXFLAGS "${CXXFLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set(FFLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "RELWITHDEBINFO")
  set(CFLAGS "${CFLAGS} ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  set(CXXFLAGS "${CXXFLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  set(FFLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "MINSIZEREL")
  set(CFLAGS "${CFLAGS} ${CMAKE_C_FLAGS_MINSIZEREL}")
  set(CXXFLAGS "${CXXFLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}")
  set(FFLAGS "${FFLAGS} ${CMAKE_Fortran_FLAGS_MINSIZEREL}")
endif()
message(STATUS "C compiler flags: " ${CFLAGS})
message(STATUS "C++ compiler flags: " ${CXXFLAGS})
message(STATUS "Fortran compiler flags: " ${FFLAGS})

# Define fine-grained targets
set_property(DIRECTORY PROPERTY EP_STEP_TARGETS
             download patch configure build install test)

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

#### Boost ####################################################################
# http://www.boost.org
ExternalProject_Add(
  boost
  PREFIX boost
  URL ${PROJECT_SOURCE_DIR}/boost_1_54_0.tar.bz2
  URL_MD5 15cb8c0803064faef0c4ddf5bc5ca279
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/boost/src/boost/bootstrap.sh
                    --prefix=${CMAKE_INSTALL_PREFIX}
                    --with-libraries=system
  INSTALL_COMMAND ${PROJECT_BINARY_DIR}/boost/src/boost/bjam
                  toolset=${BOOST_TOOLSET}
                  variant=${BOOST_BUILD_TYPE}
                  cxxflags=${CXXFLAGS}
                  cflags=${CFLAGS}
                  --debug-building
                  install
  BUILD_COMMAND ""
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### Boost.MPL extension: cartesian product ###################################
# http://www.organicvectory.com
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
ExternalProject_Add(
  hdf5ser
  PREFIX hdf5ser
  URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.11.tar.gz
  URL_MD5 1a4cc04f7dbe34e072ddcf3325717504
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5ser/src/hdf5ser/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    --prefix=${CMAKE_INSTALL_PREFIX}/hdf5ser
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### HDF5 parallel #############################################################
# http://www.hdfgroup.org/HDF5
ExternalProject_Add(
  hdf5par
  PREFIX hdf5par
  URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.11.tar.gz
  URL_MD5 1a4cc04f7dbe34e072ddcf3325717504
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5par/src/hdf5par/configure
                    CC=${MPI_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    --enable-parallel
                    --prefix=${CMAKE_INSTALL_PREFIX}/hdf5par
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

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
ExternalProject_Add(
  zlib
  PREFIX zlib
  URL ${PROJECT_SOURCE_DIR}/zlib-1.2.8.tar.gz
  URL_MD5 44d667c142d7cda120332623eab69f40
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/zlib/src/zlib/configure
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_COMMAND make CC=${UNDERLYING_C_COMPILER} CFLAGS=${CFLAGS} install
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
  DEPENDS hdf5ser
  PREFIX netcdf
  URL ${PROJECT_SOURCE_DIR}/netcdf-4.3.1.1.tar.gz
  URL_MD5 275c3b839088674c2f00fb3ac264bf11
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/netcdf/src/netcdf/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    CPPFLAGS=${CXXFLAGS}
                    CPPFLAGS=-I${CMAKE_INSTALL_PREFIX}/hdf5ser/include
                    LDFLAGS=-L${CMAKE_INSTALL_PREFIX}/hdf5ser/lib
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

#### Silo #####################################################################
# https://wci.llnl.gov/codes/silo
ExternalProject_Add(
  silo
  DEPENDS hdf5ser
  PREFIX silo
  URL ${PROJECT_SOURCE_DIR}/silo-4.9.1-bsd.tar.gz
  URL_MD5 a66223e1c0e569eb4887df2aae0ab70f
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/silo/src/silo/configure
                    CC=${MPI_C_COMPILER}
                    CFLAGS=${CFLAGS}
                    --enable-optimization
                    --enable-shared
                    --disable-silex
                    --with-hdf5=${CMAKE_INSTALL_PREFIX}/hdf5ser/include,${CMAKE_INSTALL_PREFIX}/hdf5ser/lib
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### minixml ###################################################################
# http://www.minixml.org/software.php
ExternalProject_Add(
  minixml
  PREFIX minixml
  URL ${PROJECT_SOURCE_DIR}/mxml-2.8.tar.gz
  URL_MD5 d85ee6d30de053581242c4a86e79a5d2
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/minixml/src/minixml/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CXX=${UNDERLYING_CXX_COMPILER}
                    CFLAGS=${CFLAGS}
                    CXXFLAGS=${CXXFLAGS}
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### bzip2 #####################################################################
# http://www.bzip.org
ExternalProject_Add(
  bzip2
  PREFIX bzip2
  URL ${PROJECT_SOURCE_DIR}/bzip2-1.0.6.tar.gz
  URL_MD5 00b516f4704d4a7cb50a1d97e6e8e15b
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
                CC=${UNDERLYING_C_COMPILER}
                CFLAGS=${CFLAGS}
                PREFIX=${CMAKE_INSTALL_PREFIX}
                install
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### szip ######################################################################
# http://www.hdfgroup.org/ftp/lib-external/szip
ExternalProject_Add(
  szip
  PREFIX szip
  URL ${PROJECT_SOURCE_DIR}/szip-2.1.tar.gz
  URL_MD5 902f831bcefb69c6b635374424acbead
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/szip/src/szip/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CXX=${UNDERLYING_CXX_COMPILER}
                    FC=${UNDERLYING_Fortran_COMPILER}
                    CFLAGS=${CFLAGS}
                    CXXFLAGS=${CXXFLAGS}
                    FFLAGS=${FFLAGS}
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### ISOBAR ####################################################################
# http://dx.doi.org/10.1145/2287076.2287086
ExternalProject_Add(
  isobar
  PREFIX isobar
  URL ${PROJECT_SOURCE_DIR}/isobar.0.3.0.tgz
  URL_MD5 86d123b7028ed238c0e6e49c22e7272c
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/lib
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include
  BUILD_COMMAND make install
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
                  ${PROJECT_BINARY_DIR}/isobar/src/isobar/lib/libisobar.a
                  ${CMAKE_INSTALL_PREFIX}/lib
          COMMAND ${CMAKE_COMMAND} -E copy
                  ${PROJECT_BINARY_DIR}/isobar/src/isobar/include/isobar.h
                  ${CMAKE_INSTALL_PREFIX}/include
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### APLOD #####################################################################
# http://conferences.computer.org/sc/2012/papers/1000a088.pdf
ExternalProject_Add(
  aplod
  PREFIX aplod
  URL ${PROJECT_SOURCE_DIR}/aplod.0.1.tar.gz
  URL_MD5 3e9749efb42a3d4d3165e884456812bc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make CXX=${UNDERLYING_CXX_COMPILER} CXXLAGS=${CXXFLAGS}
  INSTALL_COMMAND make install prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### ADIOS #####################################################################
# https://www.olcf.ornl.gov/center-projects/adios
ExternalProject_Add(
  adios
  DEPENDS minixml zlib bzip2 szip isobar aplod hdf5ser hdf5par
  PREFIX adios
  URL ${PROJECT_SOURCE_DIR}/adios-1.6.0.tar.gz
  URL_MD5 70073584c620b92d7842250ce778847f
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/adios/src/adios/configure
                    CC=${MPI_C_COMPILER}
                    CXX=${MPI_CXX_COMPILER}
                    CFLAGS=${CFLAGS}
                    CXXLAGS=${CXXFLAGS}
                    --disable-fortran
                    --enable-research-transports
                    --enable-skel-timing
                    --with-mxml=${CMAKE_INSTALL_PREFIX}
                    --with-phdf5=${CMAKE_INSTALL_PREFIX}/hdf5par
                    --with-hdf5=${CMAKE_INSTALL_PREFIX}/hdf5ser
                    --with-zlib=${CMAKE_INSTALL_PREFIX}
                    --with-bzip2=${CMAKE_INSTALL_PREFIX}
                    --with-szip=${CMAKE_INSTALL_PREFIX}
                    --with-isobar=${CMAKE_INSTALL_PREFIX}
                    --with-aplod=${CMAKE_INSTALL_PREFIX}
                    --prefix=${CMAKE_INSTALL_PREFIX}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### RNGSSELib2 ################################################################
# http://dx.doi.org/10.1016/j.cpc.2013.04.007
# No need for MPI, use underlying compilers
ExternalProject_Add(
  rngsse2
  PREFIX rngsse2
  URL ${PROJECT_SOURCE_DIR}/rngsselib_2.tar.gz
  URL_MD5 b78981ae65c3d801e2b878f2a38084b8
  CMAKE_ARGS -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${UNDERLYING_C_COMPILER}
             -D CMAKE_C_FLAGS:STRING=${CFLAGS}
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
                    F77=${UNDERLYING_FORTRAN_COMPILER}
                    CFLAGS=-Wno-error=return-type
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
#                    F77=${UNDERLYING_FORTRAN_COMPILER}
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
