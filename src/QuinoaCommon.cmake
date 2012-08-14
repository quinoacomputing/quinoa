# Common headers/libraries for all code

set(CMAKE_VERBOSE_MAKEFILE 1)

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

include_directories($ENV{MKLROOT}/include
)

### Libraries ###

find_library(PTHREAD_LIBRARY
             NAMES pthread
             PATHS /usr/lib
)

find_library(MATH_LIBRARY
             NAMES libm.a
             PATHS /usr/lib64
)

find_library(MKL_INTEL_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(MKL_INTEL_THREAD_LIBRARY
             NAMES mkl_intel_thread
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS $ENV{MKLROOT}/lib/intel64
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
