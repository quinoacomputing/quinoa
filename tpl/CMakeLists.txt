# Build quinoa third-party Libraries

cmake_minimum_required(VERSION 2.8.5)

project(TPL C CXX)
include(ExternalProject)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Find MPI
find_package(MPI REQUIRED)

# Set install directory based on compiler ID
string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
set(TPL_INSTALL_DIR "${PROJECT_BINARY_DIR}/../install/${COMP}")
MESSAGE(STATUS "Third-party library install dir: " ${TPL_INSTALL_DIR})

# Trilinos
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

# HDF5
ExternalProject_Add(
    hdf5
    PREFIX hdf5
    URL ${PROJECT_SOURCE_DIR}/hdf5-1.8.11.tar.gz
    URL_MD5 1a4cc04f7dbe34e072ddcf3325717504
    CONFIGURE_COMMAND CC=${MPI_C_COMPILER}
                      ${PROJECT_BINARY_DIR}/hdf5/src/hdf5/configure
                      --prefix=${TPL_INSTALL_DIR}
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 0
    LOG_BUILD 1
    LOG_INSTALL 1
)

# H5Part
ExternalProject_Add(
    h5part
    DEPENDS hdf5
    PREFIX h5part
    URL ${PROJECT_SOURCE_DIR}/H5Part-1.6.6.tar.gz
    URL_MD5 327c63d198e38a12565b74cffdf1f9d7
    CONFIGURE_COMMAND CC=${MPI_C_COMPILER}
                      ${PROJECT_BINARY_DIR}/h5part/src/h5part/configure
                      --enable-parallel
                      --enable-tools
                      --with-hdf5=${TPL_INSTALL_DIR}
                      --prefix=${TPL_INSTALL_DIR}
    BUILD_IN_SOURCE 1
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 0
    LOG_BUILD 1
    LOG_INSTALL 1
)
