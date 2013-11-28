# Third-party libraries paths

#### MKL (optional)
message(STATUS "Check for optional MKL (Intel Math Kernel Library)")

# Add MKLROOT variable to cache, initialize with MKLROOT environment variable
set(MKLROOT $ENV{MKLROOT} CACHE STRING "Root of optional MKL library. Clear this variable to disable MKL.")

message(STATUS "\tMKLROOT = ${MKLROOT}")
set(MKL_SEARCH_PATH)
list(APPEND MKL_SEARCH_PATH ${MKLROOT}/lib/intel64)

# Attempt to find libraries
set(MKL_INTERFACE_LIBRARY "NOTFOUND")
find_library(MKL_INTERFACE_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS ${MKL_SEARCH_PATH})

set(MKL_THREAD_LIBRARY "NOTFOUND")
find_library(MKL_THREAD_LIBRARY
             NAMES mkl_intel_thread
             PATHS ${MKL_SEARCH_PATH})

set(MKL_CORE_LIBRARY "NOTFOUND")
find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS ${MKL_SEARCH_PATH})

set(INTEL_OMP_RUNTIME_LIBRARY "NOTFOUND")
find_library(INTEL_OMP_RUNTIME_LIBRARY
             NAMES iomp5
             PATHS ${MKL_SEARCH_PATH}/../../../compiler/lib/intel64)

# Echo find libraries status
if (MKL_INTERFACE_LIBRARY)
  message(STATUS "\tFound MKL interface library '${MKL_INTERFACE_LIBRARY}'")
else()
  set(MKL_INTERFACE_LIBRARY "")
  message(STATUS "\tCould not find MKL interface library 'mkl_intel_ilp64'")
endif()

if (MKL_THREAD_LIBRARY)
  message(STATUS "\tFound MKL thread library '${MKL_THREAD_LIBRARY}'")
else()
  set(MKL_THREAD_LIBRARY "")
  message(STATUS "\tCould not find MKL thread library 'mkl_intel_thread'")
endif()

if (MKL_CORE_LIBRARY)
  message(STATUS "\tFound MKL core library '${MKL_CORE_LIBRARY}'")
else()
  set(MKL_CORE_LIBRARY "")
  message(STATUS "\tCould not find MKL core library 'mkl_core'")
endif()

if (INTEL_OMP_RUNTIME_LIBRARY)
  message(STATUS "\tFound Intel OpenMP runtime library '${INTEL_OMP_RUNTIME_LIBRARY}'")
else()
  set(INTEL_OMP_RUNTIME_LIBRARY "")
  message(STATUS "\tCould not find Intel OpenMP runtime library 'iomp5' required by MKL thread library 'mkl_intel_thread'")
endif()

# Define HAS_MKL macro and echo MKL status
if (MKL_INTERFACE_LIBRARY AND
    MKL_THREAD_LIBRARY AND
    MKL_CORE_LIBRARY AND
    INTEL_OMP_RUNTIME_LIBRARY)
  message(STATUS "Check for optional MKL (Intel Math Kernel Library) -- works")
  set(HAS_MKL on)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMKL_ILP64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMKL_ILP64")
else()
  message(STATUS "Check for optional MKL (Intel Math Kernel Library) -- failed")
  set(HAS_MKL off)
endif()

#### Z
find_library(Z_LIBRARY
             NAMES z
             PATHS /usr/lib64
)

#### Silo
find_library(SILO_LIBRARY
             NAMES siloh5
             PATHS ${TPL_DIR}/lib
)

#### HDF5
find_library(HDF5_LIBRARY
             NAMES hdf5
             PATHS ${TPL_DIR}/lib
)

#### TestU01
find_library(TESTU01_LIBRARY
             NAMES testu01
             PATHS ${TPL_DIR}/lib
)
find_library(TESTU01_PROBDIST_LIBRARY
             NAMES probdist
             PATHS ${TPL_DIR}/lib
)
#find_library(TESTU01_MYLIB_LIBRARY
#             NAMES mylib
#             PATHS ${TPL_DIR}/lib
#)
