# Third-party libraries paths

#### MKL (optional)
message(STATUS "Check for optional MKL (Intel Math Kernel Library)")
message(STATUS "\tBased on shell environment or cmake cache variable 'MKLROOT'")

if (DEFINED ENV{MKLROOT} OR DEFINED MKLROOT)

  # Add MKLROOT library paths to the search
  set(MKL_SEARCH_PATHS)
  if (DEFINED ENV{MKLROOT})
    message(STATUS "\tFound MKLROOT env var = $ENV{MKLROOT}")
    list(APPEND MKL_SEARCH_PATHS $ENV{MKLROOT}/lib/intel64)
  endif()
  if (MKLROOT)
    message(STATUS "\tFound MKLROOT cmake var = ${MKLROOT}")
    list(APPEND MKL_SEARCH_PATHS ${MKLROOT}/lib/intel64)
  endif()
  list(APPEND MKL_SEARCH_PATHS  # other possible install locations
              ${THIRD_PARTY_PREFIX}/lib
              ${THIRD_PARTY_LIB_DIR}
              /usr/lib64
              /usr/lib)

  # Attempt to find libraries
  find_library(MKL_INTERFACE_LIBRARY
               NAMES mkl_intel_ilp64
               PATHS ${MKL_SEARCH_PATHS})
  find_library(MKL_THREAD_LIBRARY
               NAMES mkl_intel_thread
               PATHS ${MKL_SEARCH_PATHS})
  find_library(MKL_CORE_LIBRARY
               NAMES mkl_core
               PATHS ${MKL_SEARCH_PATHS})
  find_library(INTEL_OMP_RUNTIME_LIBRARY
               NAMES iomp5
               PATHS $ENV{MKLROOT}/../compiler/lib/intel64)

  # Echo find libraries status
  if (MKL_INTERFACE_LIBRARY)
    message(STATUS "\tFound MKL interface library '${MKL_INTERFACE_LIBRARY}'")
  else()
    message(STATUS "\tCould NOT find MKL interface library 'mkl_intel_ilp64'")
  endif()

  if (MKL_THREAD_LIBRARY)
    message(STATUS "\tFound MKL thread library '${MKL_THREAD_LIBRARY}'")
  else()
    message(STATUS "\tCould NOT find MKL thread library 'mkl_intel_thread'")
  endif()

  if (MKL_CORE_LIBRARY)
    message(STATUS "\tFound MKL core library '${MKL_CORE_LIBRARY}'")
  else()
    message(STATUS "\tCould NOT find MKL core library 'mkl_core'")
  endif()

  if (INTEL_OMP_RUNTIME_LIBRARY)
    message(STATUS "\tFound Intel OpenMP runtime library '${INTEL_OMP_RUNTIME_LIBRARY}'")
  else()
    message(STATUS "\tCould NOT find Intel OpenMP runtime library 'iomp5' required by MKL thread library 'mkl_intel_thread'")
  endif()

  # Define HAS_MKL macro and echo final MKL status
  if (MKL_INTERFACE_LIBRARY AND
      MKL_THREAD_LIBRARY AND
      MKL_CORE_LIBRARY AND
      INTEL_OMP_RUNTIME_LIBRARY)
    message(STATUS "Check for optional MKL (Intel Math Kernel Library) -- works")
    set(HAS_MKL on)
  else()
    message(STATUS "Check for optional MKL (Intel Math Kernel Library) -- failed")
    set(HAS_MKL off)
  endif()
endif()


#### Pthread
find_library(PTHREAD_LIBRARY
             NAMES pthread
             PATHS /usr/lib
             PATHS /usr/lib/x86_64-redhat-linux5E/lib64
)

#### Math
find_library(MATH_LIBRARY
             NAMES m
             PATHS /usr/lib64
             PATHS /usr/lib/x86_64-redhat-linux5E/lib64
)

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
#find_library(TESTU01_PROBDIST_LIBRARY
#             NAMES probdist
#             PATHS ${TPL_DIR}/lib
#)
#find_library(TESTU01_MYLIB_LIBRARY
#             NAMES mylib
#             PATHS ${TPL_DIR}/lib
#)
