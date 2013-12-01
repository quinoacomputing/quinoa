//******************************************************************************
/*!
  \file      src/RNG/MKLRNGWrappers.C
  \author    J. Bakosi
  \date      Sat 30 Nov 2013 10:51:36 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL RNG wrappers
  \details   MKL RNG wrappers
*/
//******************************************************************************

#include <memory>
#include <vector>

#include <MKLRNGWrappers.h>

namespace rngtest {

std::vector< std::unique_ptr<tk::RNG> > g_rng;       //!< RNGs
Rsize g_rid;                                         //!< RNG id
std::vector< int > g_tid;                            //!< Global thread ids

double MKLRNGUniform(void*, void*)
//******************************************************************************
//  TestU01 MKL uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
{
  double r;
  g_rng[ g_rid ]->uniform( g_tid[g_rid], 1, &r );
  return r;
}

} // rngtest::
