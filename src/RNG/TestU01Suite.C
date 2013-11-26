//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.C
  \author    J. Bakosi
  \date      Mon 25 Nov 2013 10:55:30 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

extern "C" {
  #include <smarsa.h>
}

#include <TestU01Suite.h>

using rngtest::TestU01Suite;

double
TestU01Suite::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res )
//******************************************************************************
//  Run BirthdaySpacings test
//! \author  J. Bakosi
//******************************************************************************
{
#ifdef USE_LONGLONG
  smarsa_BirthdaySpacings( gen, res, 1, 5 * MILLION, 0, 1073741824, 2, 1 );
#else
  smarsa_BirthdaySpacings( gen, res, 10, MILLION / 2, 0, 67108864, 2, 1 );
#endif

  return res->pVal2;
}
