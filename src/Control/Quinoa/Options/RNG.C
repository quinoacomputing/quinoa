//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/RNG.C
  \author    J. Bakosi
  \date      Mon 04 Nov 2013 10:42:22 PM MST
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

unsigned int
RNG::seed( RNGType rng, const std::map< RNGType, unsigned int >& seed ) const
//******************************************************************************
//  Return seed value for RNG
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = seed.find( rng );

  if ( it != seed.end() ) {
    return it->second;          // user has specified it
  } else {
    return 0;                   // user has not specified it
  }
}

void
RNG::initFactory( RNGFactory& f,
                  std::list< std::string >& reg,
                  int nthreads,
                  const std::map< RNGType, unsigned int >& seedmap ) const
//******************************************************************************
//  Register random number generators into factory
//! \author  J. Bakosi
//******************************************************************************
{
  //! Lambda to register a random number generator into factory
  auto regRNG = [&]( RNGType rng ) {
    reg.push_back(add<MKLRNG>(f, rng, nthreads, param(rng), seed(rng,seedmap)));
  };

  regRNG( RNGType::MKL_MCG31 );
  regRNG( RNGType::MKL_R250 );
  regRNG( RNGType::MKL_MRG32K3A );
  regRNG( RNGType::MKL_MCG59 );
  regRNG( RNGType::MKL_WH );
  regRNG( RNGType::MKL_MT19937 );
  regRNG( RNGType::MKL_MT2203 );
  regRNG( RNGType::MKL_SFMT19937 );
  regRNG( RNGType::MKL_SOBOL );
  regRNG( RNGType::MKL_NIEDERR );
  regRNG( RNGType::MKL_IABSTRACT );
  regRNG( RNGType::MKL_DABSTRACT );
  regRNG( RNGType::MKL_SABSTRACT );
  regRNG( RNGType::MKL_NONDETERM );
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
