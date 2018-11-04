################################################################################
#
# \file      cmake/TPLs.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Find the third-party libraries required to build Quinoa
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
find_package(Charm)

#### MKL (optional)
find_package(MKL)
if(MKL_FOUND)
  set(HAS_MKL true)  # will become compiler define in Main/QuinoaConfig.h
endif()

#### BLAS/LAPACK library with LAPACKE C-interface
if (NOT MKL_FOUND)    # Prefer Intel's MKL for BLAS/LAPACK if available
  find_package(LAPACKE)
endif()

#### Boost
set(BOOST_INCLUDEDIR ${TPL_DIR}/include) # prefer ours
find_package(Boost 1.56.0)
if(Boost_FOUND)
  message(STATUS "Boost at ${Boost_INCLUDE_DIR} (include)")
  include_directories(${Boost_INCLUDE_DIR})
endif()

#### TUT
find_package(TUT)

#### PStreams
set(PSTREAMS_ROOT ${TPL_DIR}) # prefer ours
find_package(PStreams)

#### Hypre
find_package(Hypre 2.9.0)

#### PugiXML
set(PUGIXML_ROOT ${TPL_DIR}) # prefer ours
find_package(Pugixml)

#### PEGTL
find_package(PEGTL 2.0.0)

#### Random123
find_package(Random123)

#### RNGSSE2 library
if(ARCH MATCHES "x86")
  find_package(RNGSSE2)
endif()
if(RNGSSE2_FOUND)
  set(HAS_RNGSSE2 true)  # will become compiler define in Main/QuinoaConfig.h
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

if (NOT HDF5_FOUND)
  set(HDF5_INCLUDE_DIRS "")
endif()

#### H5Part
find_package(H5Part)

#### AEC (only for static link)
if(NOT BUILD_SHARED_LIBS)
  find_package(AEC)
endif()

#### Zlib (only for static link)
if(NOT BUILD_SHARED_LIBS AND NOT ARCH MATCHES "ppc64")
  find_package(ZLIB)
endif()

#### Zoltan2 library
find_package(Zoltan2)

#### NumDiff executable
find_package(NumDiff)

#### ExodusII library
find_package(SEACASExodus)
set(EXODUS_ROOT ${TPL_DIR}) # prefer ours
find_package(Exodiff)

#### TestU01 library
set(TESTU01_ROOT ${TPL_DIR}) # prefer ours
find_package(TestU01)
if(TestU01_FOUND)
  set(HAS_TESTU01 true)  # will become compiler define in Main/QuinoaConfig.h
endif()

### Root library
find_package(Root COMPONENTS RIO Core Tree Hist)
if (Root_FOUND)
  set(HAS_ROOT true)  # will become compiler define in Main/QuinoaConfig.h
  # Root does not support libc++ on linux, so remove if configured
  if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    string(FIND "${CMAKE_CXX_FLAGS}" "-stdlib=libc++" pos)
    if (NOT "${pos}" STREQUAL "-1")
      message(STATUS "Removing C++ compiler flag '-stdlib=libc++' as Root does not support it")
      string(REPLACE "-stdlib=libc++" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    endif()
  endif()
endif()

#### Configure Backward-cpp
set(BACKWARD_ROOT ${TPL_DIR}) # prefer ours
find_package(BackwardCpp)
if(BACKWARDCPP_FOUND)
  list(APPEND CMAKE_MODULE_PATH "${BACKWARDCPP_CMAKE_CONFIG_DIRS}")
  include(BackwardConfig)
  set(HAS_BACKWARD true)  # will become compiler define in Main/QuinoaConfig.h
  message(STATUS "Backward-cpp config: ${BACKWARD_DEFINITIONS}")
  message(STATUS "Backward-cpp libraries: ${BACKWARD_LIBRARIES}")
else()
  set(BACKWARDCPP_INCLUDE_DIRS "")
  set(BACKWARD_LIBRARIES "")
endif()

#### Configure Omega_h
find_package(Omega_h)
if(OMEGA_H_FOUND)
  set(HAS_OMEGA_H true)  # will become compiler define in Main/QuinoaConfig.h
else()
  set(OMEGA_H_INCLUDE_DIRS "")
  set(OMEGA_H_LIBRARIES "")
endif()

#### Configure HighwayHash
set(HIGHWAYHASH_ROOT ${TPL_DIR}) # prefer ours
find_package(HighwayHash)

#### Configure Brigand
set(BRIGAND_ROOT ${TPL_DIR}) # prefer ours
find_package(Brigand)

message(STATUS "------------------------------------------")

# Enable individual executables based on required TPLs found

if (CHARM_FOUND AND PUGIXML_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND
    HDF5_FOUND AND BRIGAND_FOUND AND TUT_FOUND AND PEGTL_FOUND AND Boost_FOUND)
  set(ENABLE_UNITTEST "true")
  set(UNITTEST_EXECUTABLE unittest)
else()
  message(STATUS "Executable 'unittest' will NOT be configured")
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND HYPRE_FOUND AND
    Zoltan2_FOUND AND HDF5_FOUND AND HDF5_FOUND AND BRIGAND_FOUND AND
    PEGTL_FOUND AND (MKL_FOUND OR LAPACKE_FOUND) AND Boost_FOUND)
  set(ENABLE_INCITER "true")
  set(INCITER_EXECUTABLE inciter)
else()
  message(STATUS "Executable 'inciter' will NOT be configured")
endif()

if (CHARM_FOUND AND TESTU01_FOUND AND BRIGAND_FOUND AND PEGTL_FOUND AND
    RANDOM123_FOUND AND Boost_FOUND)
  set(ENABLE_RNGTEST "true")
  set(RNGTEST_EXECUTABLE rngtest)
  set(RNGTEST_SRC_DIR ${QUINOA_SOURCE_DIR}/RNGTest)
  set(RNGTEST_BIN_DIR ${PROJECT_BINARY_DIR}/RNGTest)
else()
  message(STATUS "Executable 'rngtest' will NOT be configured")
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND PEGTL_FOUND AND
    PUGIXML_FOUND AND HDF5_FOUND AND Boost_FOUND)
  set(ENABLE_MESHCONV "true")
  set(MESHCONV_EXECUTABLE meshconv)
else()
  message(STATUS "Executable 'meshconv' will NOT be configured")
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND PEGTL_FOUND AND
    BRIGAND_FOUND AND HDF5_FOUND AND RANDOM123_FOUND AND Boost_FOUND AND
    (MKL_FOUND OR LAPACKE_FOUND))
  set(ENABLE_WALKER "true")
  set(WALKER_EXECUTABLE walker)
else()
  message(STATUS "Executable 'walker' will NOT be configured")
endif()

if (CHARM_FOUND AND SEACASExodus_FOUND AND EXODIFF_FOUND AND ROOT_FOUND
    AND PEGTL_FOUND AND PUGIXML_FOUND AND HDF5_FOUND AND Boost_FOUND)
  set(ENABLE_FILECONV "true")
  set(FILECONV_EXECUTABLE fileconv)
else()
  message(STATUS "Executable 'fileconv' will NOT be configured")
endif()
