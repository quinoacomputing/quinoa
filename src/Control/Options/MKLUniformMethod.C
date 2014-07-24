//******************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.C
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 09:19:02 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Options/MKLUniformMethod.h>

using tk::ctr::MKLUniformMethod;

const MKLUniformMethod::ParamType&
MKLUniformMethod::param( MKLUniformMethodType m ) const
//******************************************************************************
//  Return parameter based on enum
//! \author  J. Bakosi
//******************************************************************************
{
  auto it = method.find( m );

  Assert( it != method.end(),
          std::string("Cannot find parameter for MKLUniformMethod \"")
          << m << "\"" );

  return it->second;
}
