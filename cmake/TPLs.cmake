################################################################################
#
# \file      cmake/TPLs.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the third-party libraries required to build Quinoa
# \date      Fri 20 Jan 2017 01:31:49 PM MST
#
################################################################################

# Add TPL_DIR to modules directory for TPLs that provide cmake FIND_PACKAGE
# code, such as Trilinos
SET(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

# Include support for multiarch path names
include(GNUInstallDirs)

#### TPLs we attempt to find on the system #####################################

message(STATUS "------------------------------------------")

#### Charm++
set(CHARM_ROOT ${TPL_DIR}/charm)
find_package(Charm REQUIRED)

#### MKL (optional)
find_package(MKL)
if(MKL_FOUND)
  set(HAS_MKL true)  # will become compiler define in Main/QuinoaConfig.h
endif()

#### BLAS/LAPACK library with LAPACKE C-interface
if (NOT MKL_FOUND)    # Prefer Intel's MKL for BLAS/LAPACK if available
  find_package(LAPACKE REQUIRED)
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

set(CARTESIAN_PRODUCT_ROOT ${TPL_DIR}) # prefer ours
find_package(CartesianProduct REQUIRED)

#### TUT
set(TUT_ROOT ${TPL_DIR}) # prefer ours
find_package(TUT REQUIRED)

#### PStreams
set(PSTREAMS_ROOT ${TPL_DIR}) # prefer ours
find_package(PStreams REQUIRED)

#### Hypre
set(HYPRE_ROOT ${TPL_DIR}) # prefer ours
find_package(Hypre REQUIRED)

#### PugiXML
set(PUGIXML_ROOT ${TPL_DIR}) # prefer ours
find_package(Pugixml REQUIRED)

#### PEGTL
set(PEGTL_ROOT ${TPL_DIR}) # prefer ours
find_package(PEGTL REQUIRED)

#### Random123
set(Random123_ROOT ${TPL_DIR}) # prefer ours
find_package(Random123 REQUIRED)

#### RNGSSE2 library
set(RNGSSE2_ROOT ${TPL_DIR}) # prefer ours
find_package(RNGSSE2)
if(RNGSSE2_FOUND AND NOT NO_RNGSSE2)
  set(HAS_RNGSSE2 true)  # will become compiler define in Main/QuinoaConfig.h
endif()

# Error out if not a single RNG library has been found
if (NOT MKL_FOUND AND NOT Random123_FOUND AND NOT RNGSSE2_FOUND)
  message(FATAL "At least one of MKL, RNGSSE2, Random123 is required.")
endif()

### HDF5/NetCDF (NetCDF only for static link)
if(NOT BUILD_SHARED_LIBS)
  set(HDF5_PREFER_PARALLEL true)
  set(HDF5_USE_STATIC_LIBRARIES true)
  find_package(HDF5 COMPONENTS C HL)
  find_package(NetCDF)
else()
  set(HDF5_PREFER_PARALLEL true)
  find_package(HDF5 COMPONENTS C HL)
endif()

#### H5Part
set(H5PART_ROOT ${TPL_DIR}) # prefer ours
find_package(H5Part REQUIRED)

#### AEC (only for static link)
if(NOT BUILD_SHARED_LIBS)
  set(AEC_ROOT ${TPL_DIR}) # prefer ours
  find_package(AEC REQUIRED)
endif()

#### Zlib (only for static link)
if(NOT BUILD_SHARED_LIBS)
  find_package(ZLIB REQUIRED)
endif()

#### Zoltan2 library
find_package(Zoltan2 REQUIRED)

#### ExodusII library
find_package(SEACASExodus REQUIRED)
set(EXODUS_ROOT ${TPL_DIR}) # prefer ours
find_package(Exodiff REQUIRED)

#### TestU01 library
set(TESTU01_ROOT ${TPL_DIR}) # prefer ours
find_package(TestU01)
if(TestU01_FOUND AND NOT NO_TESTU01)
  set(HAS_TESTU01 true)  # will become compiler define in Main/QuinoaConfig.h
endif()

message(STATUS "------------------------------------------")
