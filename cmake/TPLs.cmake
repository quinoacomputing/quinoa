# Find third-party libraries

# Add TPL_DIR/include to modules directory for TPLs that provide cmake
# FIND_PACKAGE code, such as Trilinos
SET(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

#### TPLs we attempt to find on the system #####################################

#### MKL (optional)
find_package(MKL)
if(MKL_FOUND)
  set(HAS_MKL true)  # will become compiler define in Main/Config.h
endif()

#### BLAS/LAPACK library with LAPACKE C-interface
if (NOT MKL_FOUND)    # Prefer Intel's MKL for BLAS/LAPACK if available
  # If MKL is unavailable, prefer ours then fall back to system
  find_path(LAPACKE_PATH lapacke.h DOC "C-interface to LAPACK")
  find_library(LAPACKE_LIBRARY NAMES lapacke HINTS ${TPL_DIR}/lib NO_DEFAULT_PATH)
  find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
  message(STATUS "Found LAPACK C-interface ${LAPACKE_PATH}/lapacke.h and ${LAPACKE_LIBRARY}")
endif()

#### Boost
if(NOT NO_SYSTEM_BOOST)
  set(BOOST_INCLUDEDIR ${TPL_DIR}/include) # prefer ours
  find_package(Boost REQUIRED)
endif()
if(Boost_FOUND)
  message(STATUS "Boost at ${Boost_INCLUDE_DIR} (include)")
  include_directories(${Boost_INCLUDE_DIR})
endif()

#### PStreams
set(PSTREAMS_ROOT ${TPL_DIR}) # prefer ours
find_package(PStreams REQUIRED)

#### Hypre
set(HYPRE_ROOT ${TPL_DIR}) # prefer ours
find_package(Hypre REQUIRED)

#### PugiXML
set(PUGIXML_ROOT ${TPL_DIR}) # prefer ours
find_package(pugixml REQUIRED)

#### HDF5 (only for static link)
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_PREFER_PARALLEL true)
  set(HDF5_USE_STATIC_LIBRARIES true)
  find_package(HDF5 COMPONENTS HL)
  find_package(NetCDF)
endif()

#### AEC (only for static link)
if(NOT BUILD_SHARED_LIBS)
  find_package(AEC REQUIRED)
endif()

#### TPLs we always want ours ##################################################

#### Zoltan2 library
find_package(Zoltan2 REQUIRED)
if(Zoltan2_FOUND)
  message(STATUS "Found Zoltan2: ${Zoltan2_LIBRARY_DIRS}")
endif()

#### ExodusII library
find_package(SEACASExodus REQUIRED)
if(SEACASExodus_FOUND)
  message(STATUS "Found SEACASExodus: ${SEACASExodus_LIBRARY_DIRS}")
endif()
find_package(SEACASNemesis REQUIRED)
if(SEACASNemesis_FOUND)
  message(STATUS "Found SEACASNemesis: ${SEACASNemesis_LIBRARY_DIRS}")
endif()
find_package(SEACASExodiff REQUIRED)
if(SEACASExodiff_FOUND)
  message(STATUS "Found SEACASExodiff: ${SEACASExodiff_LIBRARY_DIRS}")
endif()

#### RNGSSE2 library
set(RNGSSE_LIBRARY "NOTFOUND")
find_library(RNGSSE_LIBRARY
             NAMES rngsse
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)

#### TestU01 library
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
set(TESTU01_MYLIB_LIBRARY "NOTFOUND")
find_library(TESTU01_MYLIB_LIBRARY
             NAMES mylib
             PATHS ${TPL_DIR}/lib
             NO_DEFAULT_PATH
             REQUIRED)
