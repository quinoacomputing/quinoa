# Headers and libraries paths

#### MKL
find_library(MKL_INTERFACE_LIBRARY
             NAMES mkl_intel_ilp64
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(MKL_THREAD_LIBRARY
             NAMES mkl_intel_thread
             PATHS $ENV{MKLROOT}/lib/intel64
)
find_library(MKL_CORE_LIBRARY
             NAMES mkl_core
             PATHS $ENV{MKLROOT}/lib/intel64
)

find_library(INTEL_OMP_LIBRARY
             NAMES iomp5
             PATHS $ENV{MKLROOT}/../compiler/lib/intel64
)

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
