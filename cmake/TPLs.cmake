# Third-party libraries paths

#### MKL (optional)
message(STATUS "Check for optional MKL (Intel Math Kernel Library)")

# Attempt to find MKL libraries
find_library(MKL_INTERFACE_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
             NAMES mkl_sequential
             PATHS $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_THREADED_LAYER_LIBRARY
             NAMES mkl_intel_thread
             PATHS $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS $ENV{MKLROOT}/lib/intel64
                   $ENV{INTEL}/mkl/lib/intel64
             NO_DEFAULT_PATH)

find_path(MKL_INCLUDE_PATH mkl.h
          $ENV{MKLROOT}/include
          $ENV{INTEL}/mkl/include
          NO_DEFAULT_PATH)

# Define HAS_MKL macro and echo MKL status
if (MKL_INTERFACE_LIBRARY AND
    MKL_SEQUENTIAL_LAYER_LIBRARY AND
    MKL_THREADED_LAYER_LIBRARY AND
    MKL_CORE_LIBRARY)
  message(STATUS "Check for optional MKL (Intel Math Kernel Library) -- works")
  set(HAS_MKL on)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64 -m64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64 -m64")
else()
  message(WARNING " Check for optional MKL (Intel Math Kernel Library) -- failed:\n Intel MKL VSL RNGs will not be available!")
  set(HAS_MKL off)
endif()

##### TBB (optional)
#message(STATUS "Check for optional TBB (Threading Building Blocks Library)")
#
## Add TBBROOT variable to cache, initialize with TBBROOT environment variable
#set(TBBROOT $ENV{TBBROOT} CACHE STRING "Root of optional TBB library. Clear this variable to disable TBB.")
#
#message(STATUS "  TBBROOT = ${TBBROOT}")
#set(TBB_SEARCH_PATH)
#list(APPEND TBB_SEARCH_PATH ${TBBROOT}/lib/intel64/gcc4.4)
#
## Attempt to find libraries
#set(TBB_LIBRARY "NOTFOUND")
#if (CMAKE_BUILD_TYPE MATCHES DEBUG OR CMAKE_BUILD_TYPE MATCHES RELWITHDEBINFO)
#  find_library(TBB_LIBRARY
#               NAMES tbb_debug
#               PATHS ${TBB_SEARCH_PATH}
#               NO_DEFAULT_PATH)
#else()
#  find_library(TBB_LIBRARY
#               NAMES tbb
#               PATHS ${TBB_SEARCH_PATH}
#               NO_DEFAULT_PATH)
#endif()
#
## Echo find libraries status
#if (TBB_LIBRARY)
#  message(STATUS "  Found TBB library '${TBB_LIBRARY}'")
#else()
#  set(TBB_LIBRARY "")
#  if (CMAKE_BUILD_TYPE MATCHES DEBUG OR CMAKE_BUILD_TYPE MATCHES RELWITHDEBINFO)
#    message(STATUS "  Could not find TBB library 'tbb_debug'")
#  else()
#    message(STATUS "  Could not find TBB library 'tbb'")
#  endif()
#endif()
#
## Define HAS_TBB macro and echo TBB status
#if (TBB_LIBRARY)
#  message(STATUS "Check for optional TBB (Threading Building Blocks Library) -- works")
#  set(HAS_TBB on)
#else()
#  message(STATUS "Check for optional TBB (Threading Building Blocks Library) -- failed")
#  set(HAS_TBB off)
#endif()

#### Z
#set(Z_LIBRARY "NOTFOUND")
#find_library(Z_LIBRARY
#             NAMES z
#             PATHS ${TPL_DIR}/lib
#             REQUIRED
#             NO_DEFAULT_PATH)

##### Silo
#set(SILO_LIBRARY "NOTFOUND")
#find_library(SILO_LIBRARY
#             NAMES siloh5
#             PATHS ${TPL_DIR}/lib
#             NO_DEFAULT_PATH
#             REQUIRED)

#### ExodusII
set(EXODUS_LIBRARY "NOTFOUND")
find_library(EXODUS_LIBRARY
             NAMES exodus
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### Nemesis
set(NEMESIS_LIBRARY "NOTFOUND")
find_library(NEMESIS_LIBRARY
             NAMES nemesis
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### NetCDF
set(NETCDF_LIBRARY "NOTFOUND")
find_library(NETCDF_LIBRARY
             NAMES netcdf
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### HDF5
set(HDF5_LIBRARY "NOTFOUND")
find_library(HDF5_LIBRARY
             NAMES hdf5
             PATHS ${TPL_DIR}/hdf5ser/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### RNGSSE2
set(RNGSSE_LIBRARY "NOTFOUND")
find_library(RNGSSE_LIBRARY
             NAMES rngsse
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### Boost C++ libraries
set(BOOST_INCLUDEDIR ${TPL_DIR}/include) # prefer ours
find_package(Boost)
if(Boost_FOUND)
  message(STATUS "Boost include dir ${Boost_INCLUDE_DIR}")
  include_directories(${Boost_INCLUDE_DIR})
endif()

#### TestU01
set(TESTU01_LIBRARY "NOTFOUND")
find_library(TESTU01_LIBRARY
             NAMES testu01
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)
set(TESTU01_PROBDIST_LIBRARY "NOTFOUND")
find_library(TESTU01_PROBDIST_LIBRARY
             NAMES probdist
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)
#set(TESTU01_MYLIB_LIBRARY "NOTFOUND")
#find_library(TESTU01_MYLIB_LIBRARY
#             NAMES mylib
#             PATHS ${TPL_DIR}/lib
#             NO_DEFAULT_PATH
#             REQUIRED)
