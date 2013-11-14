//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/RNG.C
  \author    J. Bakosi
  \date      Thu Nov 14 08:10:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Quinoa/Options/RNG.h>

using quinoa::ctr::RNG;

const RNG::ParamType&
RNG::param( RNGType rng ) const
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
RNG::mkl_seed( RNGType rng, const MKLRNGParameters& mklparam ) const
//******************************************************************************
//  Return seed value from MKLRNGParams
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = mklparam.find( rng );

  if ( it != mklparam.end() ) { // user has specified it
    return it->second.get< ctr::seed >();
  } else {                      // user has not specified it, return default
    return 0;
  }
}

quinoa::ctr::MKLUniformMethodType
RNG::mkl_uniform_method( RNGType rng, const MKLRNGParameters& mklparam ) const
//******************************************************************************
//  Return uniform method from MKLRNGParams
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = mklparam.find( rng );

  if ( it != mklparam.end() ) { // user has specified it
    return it->second.get< ctr::uniform_method >();
  } else {                      // user has not specified it, return default
    return MKLUniformMethodType::STANDARD;
  }
}

quinoa::ctr::MKLGaussianMethodType
RNG::mkl_gaussian_method( RNGType rng, const MKLRNGParameters& mklparam ) const
//******************************************************************************
//  Return Gaussian method from MKLRNGParams
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = mklparam.find( rng );

  if ( it != mklparam.end() ) { // user has specified it
    return it->second.get< ctr::gaussian_method >();
  } else {                      // user has not specified it, return default
    return MKLGaussianMethodType::BOXMULLER;
  }
}

RNG::LibType
RNG::lib( RNGType rng ) const
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
