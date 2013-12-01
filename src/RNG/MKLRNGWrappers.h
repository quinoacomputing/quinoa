//******************************************************************************
/*!
  \file      src/RNG/MKLRNGWrappers.h
  \author    J. Bakosi
  \date      Sat 30 Nov 2013 10:51:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL RNG wrappers
  \details   MKL RNG wrappers
*/
//******************************************************************************
#ifndef MKLRNGWrappers_h
#define MKLRNGWrappers_h

#include <RNG.h>

namespace rngtest {

using Rsize = std::vector< std::unique_ptr<tk::RNG> >::size_type;

double MKLRNGUniform(void*, void*);

} // rngtest::

#endif // MKLRNGWrappers_h
