# Headers and libraries paths

set(CMAKE_VERBOSE_MAKEFILE 1)

find_library(MKL_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(MKL_THREAD_LIBRARY
             NAMES mkl_intel_thread
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(INTEL_OMP_LIBRARY
             NAMES iomp5
             PATHS $ENV{MKLROOT}/../compiler/lib/intel64
)

find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS $ENV{MKLROOT}/lib/intel64
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

find_library(SILO_LIBRARY
             NAMES silo
             PATHS ${TPL_DIR}/lib
)
