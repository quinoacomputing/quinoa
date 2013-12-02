//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Mon 02 Dec 2013 06:03:01 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <InitRNGTest.h>
#include <Config.h>
#include <RNGTestDriver.h>

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Random number generator test suite
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< rngtest::RNGTestDriver, rngtest::echoTPL >
                 ( argc,
                   argv,
                   "Quinoa: Random number generator (RNG) test suite",
                   RNGTEST_EXECUTABLE);
}
