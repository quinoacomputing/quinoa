//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/RNG.C
  \author    J. Bakosi
  \date      Sat 09 Nov 2013 03:42:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Quinoa/Options/RNG.h>

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
RNG::seed( RNGType rng, const MKLRNGParam& mklparam ) const
//******************************************************************************
//  Return seed value for RNG
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = mklparam.find( rng );

  if ( it != mklparam.end() ) {
    return it->second.get<ctr::seed>();  // user has specified it
  } else {
    return 0;                            // user has not specified it
  }
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
