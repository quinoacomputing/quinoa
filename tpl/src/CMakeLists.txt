# Build quinoa thirdparty Libraries

cmake_minimum_required(VERSION 2.8.5)

project(QUINOA_TPL C CXX)
include(ExternalProject)
set(CMAKE_VERBOSE_MAKEFILE 1)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/")

# Get underlying compiler
include(GetUnderlyingCompiler)

ExternalProject_Add(
    zoltan
    PREFIX ${PROJECT_BINARY_DIR}
    URL ${PROJECT_SOURCE_DIR}/trilinos-11.2.4-Source.tar.gz
    URL_MD5 29abd4c5bd4bb39aef87e151bfe12541
    CMAKE_ARGS -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF
               -D Trilinos_ENABLE_Zoltan:BOOL=ON
               -D TPL_ENABLE_MPI:BOOL=ON
               -D Trilinos_ENABLE_Fortran:BOOL=OFF
               -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
               -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER}
               -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER}
               -D CMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/../install/${UNDERLYING_C_COMPILER}
)

# Third-Party library dependencies
#add_dependencies(${QUINOA_EXECUTABLE} zoltan)
