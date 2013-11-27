//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.C
  \author    J. Bakosi
  \date      Wed 27 Nov 2013 12:50:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

extern "C" {
  #include <smarsa.h>
  #include <sknuth.h>
}

#include <TestU01Suite.h>

using rngtest::TestU01Suite;

double
TestU01Suite::BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res )
//******************************************************************************
//  Run Marsdaglia's BirthdaySpacings test
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

double
TestU01Suite::Collision( unif01_Gen* gen, sknuth_Res2* res )
//******************************************************************************
//  Run Knuth's Collision test
//! \author  J. Bakosi
//******************************************************************************
{
  sknuth_Collision (gen, res, 1, 5 * MILLION, 0, 65536, 2);

  return res->Pois->pVal2;
}
