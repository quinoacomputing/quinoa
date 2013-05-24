# Common headers/libraries for all code

set(CMAKE_VERBOSE_MAKEFILE 1)

# Prefer static libraries
if(LINK_STRATEGY STREQUAL "STATIC")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

### Include Dirs ###

#find_path(EXODUS_INCLUDE_DIR exodusII.h
#          ${EXODUS_INCLUDE_DIR}
#          ${THIRD_PARTY_INCLUDE}
#          ${THIRD_PARTY_PREFIX}/include
#          ${INTEL_INCLUDE_GUESS}
#)
#
#find_path(NETCDF_INCLUDE_DIR netcdf.h
#          ${NETCDF_INCLUDE_DIR}
#          ${THIRD_PARTY_INCLUDE}
#          ${THIRD_PARTY_PREFIX}/include
#          ${INTEL_INCLUDE_GUESS}
#)

include_directories($ENV{INTEL}/mkl/include
)

### Libraries ###

find_library(MKL_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS $ENV{MKLROOT}/lib/intel64
             PATHS $ENV{INTEL}/mkl/lib/intel64
)

find_library(MKL_THREAD_LIBRARY
             NAMES mkl_intel_thread
             PATHS $ENV{MKLROOT}/lib/intel64
             PATHS $ENV{INTEL}/mkl/lib/intel64
)

# Linking MKL/OpenMP with clang needs explicit linking to openmp
if(PLATFORM MATCHES "clang")
  find_library(INTEL_OMP_LIBRARY
               NAMES iomp5
               PATHS $ENV{INTEL}/lib/intel64
  )
endif()

find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS $ENV{MKLROOT}/lib/intel64
             PATHS $ENV{INTEL}/mkl/lib/intel64
)

find_library(PTHREAD_LIBRARY
             NAMES pthread
             PATHS /usr/lib
             PATHS /usr/lib/x86_64-redhat-linux5E/lib64
)

find_library(MATH_LIBRARY
             NAMES m
             PATHS /usr/lib64
             PATHS /usr/lib/x86_64-redhat-linux5E/lib64
)

#find_library(EXO_C_LIBRARY
#             NAMES exoIIv2c
#             PATHS ${EXODUS_C_LIB_PATH}
#                   ${EXODUS_LIB_PATH}
#                   ${THIRD_PARTY_LIB_DIR}
#                   ${THIRD_PARTY_PREFIX}/lib
#                   ${INTEL_LIB_GUESS}
#)
#
#find_library(NETCDF_LIBRARY
#             NAMES libnetcdf.a
#             PATHS ${NETCDF_LIB_PATH}
#                   ${THIRD_PARTY_LIB_DIR}
#                   ${THIRD_PARTY_PREFIX}/lib
#                   NO_DEFAULT_PATH
#)
#find_library(HDF5HL_LIBRARY 
#             NAMES libhdf5_hl.a
#             PATHS ${HDF5_LIB_PATH}
#                   ${THIRD_PARTY_LIB_DIR}
#                   ${THIRD_PARTY_PREFIX}/lib
#                   ${INTEL_LIB_GUESS}
#)
#
#find_library(HDF5_LIBRARY 
#             NAMES libhdf5.a
#             PATHS ${HDF5_LIB_PATH}
#                   ${THIRD_PARTY_LIB_DIR}
#                   ${THIRD_PARTY_PREFIX}/lib
#                   ${INTEL_LIB_GUESS}
#)
