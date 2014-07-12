//******************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:10:43 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Options/MKLGaussianMethod.h>

using tk::ctr::MKLGaussianMethod;

const MKLGaussianMethod::ParamType&
MKLGaussianMethod::param( MKLGaussianMethodType m ) const
//******************************************************************************
//  Return parameter based on enum
//! \author  J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  auto it = method.find( m );

  Assert(it != method.end(),
         std::string("Cannot find parameter for MKLGaussianMethod \"") + m + "\"");

  return it->second;
}
