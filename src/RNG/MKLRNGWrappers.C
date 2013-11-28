//******************************************************************************
/*!
  \file      src/RNG/MKLRNGWrappers.C
  \author    J. Bakosi
  \date      Wed 27 Nov 2013 09:26:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL RNG wrappers
  \details   MKL RNG wrappers
*/
//******************************************************************************

#include <memory>

#include <MKLRNGWrappers.h>
#include <RNG.h>

namespace rngtest {

std::unique_ptr< tk::RNG > g_rng;       //!< RNG used by TestU01
int g_tid;                              //!< Global thread id used by TestU01

double MKLRNGUniform()
//******************************************************************************
//  TestU01 MKL uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r;
  g_rng->uniform( g_tid, 1, &r );
  return r;
}

} // rngtest::
