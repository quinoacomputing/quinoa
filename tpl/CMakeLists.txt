# Build quinoa third-party Libraries

cmake_minimum_required(VERSION 2.8.5)

project(TPL C CXX Fortran)
include(ExternalProject)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

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

# Set install directory based on compiler ID
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
set(TPL_INSTALL_DIR "${PROJECT_SOURCE_DIR}/install/${COMP}")
MESSAGE(STATUS "Third-party library install dir: " ${TPL_INSTALL_DIR})

#### PEGTL (Parsing Expression Grammar Template Library) ######################
# http://code.google.com/p/pegtl
# Header only
ExternalProject_Add(
  pegtl
  PREFIX pegtl
  URL ${PROJECT_SOURCE_DIR}/pegtl-0.32-patched.tar.gz
  URL_MD5 425c9717783b11f485483dcb66f7d60a
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
          ${TPL_INSTALL_DIR}/include/pegtl
)

#### Boost ####################################################################
# http://www.boost.org
# Header only
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
          ${TPL_INSTALL_DIR}/include/boost
)

#### Trilinos #################################################################
# http://trilinos.sandia.gov
ExternalProject_Add(
  trilinos
  PREFIX trilinos
  URL ${PROJECT_SOURCE_DIR}/trilinos-11.2.4-Source.tar.gz
  URL_MD5 29abd4c5bd4bb39aef87e151bfe12541
  CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
             -D Trilinos_ENABLE_Zoltan:BOOL=ON
             -D TPL_ENABLE_MPI:BOOL=ON
             -D Trilinos_ENABLE_Fortran:BOOL=OFF
             -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
             -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
             -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
             -D CMAKE_INSTALL_PREFIX:PATH=${TPL_INSTALL_DIR}
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### HDF5 #####################################################################
# http://www.hdfgroup.org/HDF5
ExternalProject_Add(
  hdf5
  PREFIX hdf5
  URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.11.tar.gz
  URL_MD5 1a4cc04f7dbe34e072ddcf3325717504
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/hdf5/src/hdf5/configure
                    CC=${MPI_C_COMPILER}
                    --prefix=${TPL_INSTALL_DIR}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### H5Part ###################################################################
# http://vis.lbl.gov/Research/H5Part
ExternalProject_Add(
  h5part
  DEPENDS hdf5
  PREFIX h5part
  URL ${PROJECT_SOURCE_DIR}/H5Part-1.6.6.tar.gz
  URL_MD5 327c63d198e38a12565b74cffdf1f9d7
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/h5part/src/h5part/configure
                    CC=${MPI_C_COMPILER}
                    --enable-parallel
                    --with-hdf5=${TPL_INSTALL_DIR}
                    --prefix=${TPL_INSTALL_DIR}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

#### Silo #####################################################################
# https://wci.llnl.gov/codes/silo
ExternalProject_Add(
  silo
  DEPENDS hdf5
  PREFIX silo
  URL ${PROJECT_SOURCE_DIR}/silo-4.9.1-bsd.tar.gz
  URL_MD5 a66223e1c0e569eb4887df2aae0ab70f
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/silo/src/silo/configure
                    CC=${MPI_C_COMPILER}
                    --enable-optimization
                    --disable-silex
                    --with-hdf5=${TPL_INSTALL_DIR}/include,${TPL_INSTALL_DIR}/lib
                    --prefix=${TPL_INSTALL_DIR}
  BUILD_IN_SOURCE 1
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
  URL ${PROJECT_SOURCE_DIR}/TestU01.zip
  URL_MD5 ff9bd5106178c1ca0665691db8e57293
  CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/testu01/src/testu01/configure
                    CC=${UNDERLYING_C_COMPILER}
                    CXX=${UNDERLYING_CXX_COMPILER}
                    F77=${UNDERLYING_FORTRAN_COMPILER}
                    CFLAGS=-Wno-error=return-type
                    --prefix=${TPL_INSTALL_DIR}
  BUILD_IN_SOURCE 1
  LOG_DOWNLOAD 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)
