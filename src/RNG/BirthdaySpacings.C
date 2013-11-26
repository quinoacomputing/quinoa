//******************************************************************************
/*!
  \file      src/RNG/BirthdaySpacings.C
  \author    J. Bakosi
  \date      Mon 25 Nov 2013 09:12:40 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical tests suggested by George Marsaglia
  \details   Statistical tests suggested by George Marsaglia
*/
//******************************************************************************

#include <iostream>

extern "C" {
  #include <smarsa.h>
}

#include <BirthdaySpacings.h>

using rngtest::BirthdaySpacings;

double
BirthdaySpacings::run()
//******************************************************************************
//  Run BirthdaySpacings test
//! \author  J. Bakosi
//******************************************************************************
{
  // Pretty awful that TestU01 does not guarantee the constness of gen
  unif01_Gen* gen = const_cast< unif01_Gen* >( m_gen );

  PoissonResPtr::element_type* res = m_res.get();

#ifdef USE_LONGLONG
  smarsa_BirthdaySpacings( gen, res, 1, 5 * MILLION, 0, 1073741824, 2, 1 );
#else
  smarsa_BirthdaySpacings( gen, res, 10, MILLION / 2, 0, 67108864, 2, 1 );
#endif

  return m_res->pVal2;
}
