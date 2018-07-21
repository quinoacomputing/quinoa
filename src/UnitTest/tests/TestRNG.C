// *****************************************************************************
/*!
  \file      src/UnitTest/tests/TestRNG.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tess for directory RNG
  \details   Unit test for directory RNG.
*/
// *****************************************************************************

#include "TUTConfig.h"

#include "tests/RNG/TestRNG.h"

#ifdef HAS_MKL
  #include "tests/RNG/TestMKLRNG.h"
#endif

#ifdef HAS_RNGSSE2
  #include "tests/RNG/TestRNGSSE.h"
#endif

#include "tests/RNG/TestRandom123.h"
