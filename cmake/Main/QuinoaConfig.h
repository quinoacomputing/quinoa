// *****************************************************************************
/*!
  \file      src/Main/QuinoaConfig.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     CMake header input file with information to be passed to cmake
  \details   CMake header input file with information to be passed to cmake. See
    also src/CMakeLists.txt.
 */
// *****************************************************************************
#ifndef QuinoaConfig_h
#define QuinoaConfig_h

#define INCITER_EXECUTABLE           "inciter"
#define RNGTEST_EXECUTABLE           "rngtest"
#define UNITTEST_EXECUTABLE          "unittest"
#define MESHCONV_EXECUTABLE          "meshconv"
#define WALKER_EXECUTABLE            "walker"

namespace tk {
 
#define VERSION                      "0.1 (C16015)"
#define GIT_COMMIT                   "Quinoa_v0.1-885-g7f160cd"
#define MPI_COMPILER                 "/usr/bin/mpicxx"
#define COMPILER                     "/usr/bin/g++"
#define BUILD_HOSTNAME               "aditya-VirtualBox"
#define BUILD_TYPE                   "DEBUG"
#define BUILD_DATE                   "Mon Jun  5 15:56:33 MDT 2017"
#define REGRESSION_DIR               "/home/aditya/quinoa/src/../regression"

// Compile-time options

// Data layout for particle data
#define PARTICLE_DATA_LAYOUT_AS_PARTICLE_MAJOR
/* #undef PARTICLE_DATA_LAYOUT_AS_EQUATION_MAJOR */

// Data layout for mesh node data
#define FIELD_DATA_LAYOUT_AS_FIELD_MAJOR
/* #undef FIELD_DATA_LAYOUT_AS_EQUATION_MAJOR */

// Optional TPLs
/* #undef HAS_MKL */
#define HAS_RNGSSE2
#define HAS_TESTU01

} // tk::

#endif // QuinoaConfig_h
