//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Wed Mar 19 13:16:52 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <RNGTestDriver.h>

namespace rngtest {

void echoTPL(const tk::Print& /*print*/)
//******************************************************************************
//  Echo TPL version informaion for libs specific to RNGTest
//! \author  J. Bakosi
//******************************************************************************
{
}

} // rngtest::

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Random number generator test suite
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< rngtest::RNGTestDriver >
                 ( argc,
                   argv,
                   "Quinoa: Random number generator (RNG) test suite",
                   RNGTEST_EXECUTABLE,
                   rngtest::echoTPL );
}
