// *****************************************************************************
/*!
  \file      src/UnitTest/tests/TestControl.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tess for directory Control
  \details   Unit test for directory Control.
*/
// *****************************************************************************

#include "TUTConfig.h"

#include "tests/Control/TestSystemComponents.h"
#include "tests/Control/TestControl.h"
#include "tests/Control/TestFileParser.h"
#include "tests/Control/TestStringParser.h"
#include "tests/Control/TestToggle.h"

#ifdef HAS_MKL
  #include "tests/Control/Options/TestMKLUniformMethod.h"
  #include "tests/Control/Options/TestMKLGaussianMethod.h"
  #include "tests/Control/Options/TestMKLBetaMethod.h"
  #include "tests/Control/Options/TestMKLGammaMethod.h"
#endif

#include "tests/Control/Options/TestRNG.h"
