//******************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:10:33 PM MDT
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
  using tk::operator+;

  auto it = method.find( m );

  Assert( it != method.end(),
          std::string("Cannot find parameter for MKLUniformMethod \"") + m + "\"");

  return it->second;
}
