//******************************************************************************
/*!
  \file      src/Control/Options/RNG.C
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 03:36:24 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Options/RNG.h>

using tk::ctr::RNG;

const RNG::ParamType&
RNG::param( RNGType rng ) const
//******************************************************************************
//  Return parameter based on enum
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = brng.find( rng );

  Assert( it != end(brng),
          std::string("Cannot find parameter for RNG \"") << rng << "\"" );

  return it->second;
}

RNG::LibType
RNG::lib( RNGType rng ) const
//******************************************************************************
//  Return RNG library type based on Enum
//! \author  J. Bakosi
//******************************************************************************
{
  const auto& n = name( rng );

  if ( found( "MKL", n ) ) return RNGLibType::MKL;
  else if ( found( "RNGSSE", n ) ) return RNGLibType::RNGSSE;
  else if ( found( "PRAND", n) ) return RNGLibType::PRAND;
  else return RNGLibType::NO_LIB;
}

bool
RNG::found( const std::string& kw, const std::string& str ) const
//******************************************************************************
//  Search for 'kw' in 'str'
//! \author  J. Bakosi
//******************************************************************************
{
  std::string::size_type f = str.find( kw );
  return f != std::string::npos ? true : false;
}
