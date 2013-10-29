//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/RNG.C
  \author    J. Bakosi
  \date      Tue Oct 29 15:17:38 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <boost/functional/factory.hpp>

#include <Quinoa/Options/RNG.h>
#include <MKLRNG.h>

using namespace quinoa::ctr;

const RNG::ParamType&
RNG::param(RNGType rng) const
//******************************************************************************
//  Return parameter based on enum
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  auto it = brng.find(rng);

  Assert(it != brng.end(), tk::ExceptType::FATAL,
         std::string("Cannot find parameter for RNG \"") + rng + "\"");

  return it->second;
}

void
RNG::initFactory( RNGFactory& f, std::list<std::string>& names,
                  int nthreads, unsigned int seed ) const
//******************************************************************************
//  Register random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
 names.push_back( add<MKLRNG>(f, RNGType::MKL_MCG31,
                              nthreads, param(RNGType::MKL_MCG31), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_R250,
                              nthreads, param(RNGType::MKL_R250), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_MRG32K3A,
                              nthreads, param(RNGType::MKL_MRG32K3A), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_MCG59,
                              nthreads, param(RNGType::MKL_MCG59), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_WH,
                              nthreads, param(RNGType::MKL_WH), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_MT19937,
                              nthreads, param(RNGType::MKL_MT19937), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_MT2203,
                              nthreads, param(RNGType::MKL_MT2203), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_SFMT19937,
                              nthreads, param(RNGType::MKL_SFMT19937), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_SOBOL,
                              nthreads, param(RNGType::MKL_SOBOL), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_NIEDERR,
                              nthreads, param(RNGType::MKL_NIEDERR), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_IABSTRACT,
                              nthreads, param(RNGType::MKL_IABSTRACT), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_DABSTRACT,
                              nthreads, param(RNGType::MKL_DABSTRACT), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_SABSTRACT,
                              nthreads, param(RNGType::MKL_SABSTRACT), seed) );

 names.push_back( add<MKLRNG>(f, RNGType::MKL_NONDETERM,
                              nthreads, param(RNGType::MKL_NONDETERM), seed) );
}

RNG::LibType
RNG::lib(RNGType rng) const
//******************************************************************************
//  Return RNG library tpe based on Enum
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  auto it = names.find(rng);

  Assert(it != names.end(), tk::ExceptType::FATAL,
         std::string("Cannot find name for RNG \"") + rng + "\"");

  if (found("MKL", it->second)) {
    return RNGLibType::MKL;
  } else if (found("RNGSSELIB", it->second)) {
    return RNGLibType::RNGSSELIB;
  } else if (found("PRAND", it->second)) {
    return RNGLibType::PRAND;
  } else {
    return RNGLibType::NO_LIB;
  }
}

bool
RNG::found(const std::string& kw, const std::string& str) const
//******************************************************************************
//  Search for 'kw' in 'str'
//! \author  J. Bakosi
//******************************************************************************
{
  std::size_t f = str.find(kw);

  if (f != std::string::npos) {
    return true;
  } else {
    return false;
  }
}
