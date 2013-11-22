//******************************************************************************
/*!
  \file      src/RNG/MKLRNGWrappers.C
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 06:15:11 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL RNG wrappers
  \details   MKL RNG wrappers
*/
//******************************************************************************

#include <memory>

#include <MKLRNGWrappers.h>
#include <RNG.h>

namespace rngtest {

std::unique_ptr< tk::RNG > g_rng;

double MKLRNGUniform()
//******************************************************************************
//  TestU01 MKL uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r;
  g_rng->uniform( 0, 1, &r );
  return r;
}

} // rngtest::
