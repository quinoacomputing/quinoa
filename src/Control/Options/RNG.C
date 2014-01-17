//******************************************************************************
/*!
  \file      src/Control/Options/RNG.C
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 08:37:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Options/RNG.h>

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

RNG::LibType
RNG::lib( RNGType rng ) const
//******************************************************************************
//  Return RNG library type based on Enum
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  auto it = names.find(rng);

  Assert(it != names.end(), tk::ExceptType::FATAL,
         std::string("Cannot find name for RNG \"") + rng + "\"");

  if (found("MKL", it->second)) {
    return RNGLibType::MKL;
  } else if (found("RNGSSE", it->second)) {
    return RNGLibType::RNGSSE;
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
  std::string::size_type f = str.find(kw);

  if (f != std::string::npos) {
    return true;
  } else {
    return false;
  }
}
