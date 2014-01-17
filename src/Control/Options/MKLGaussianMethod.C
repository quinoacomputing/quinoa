//******************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.C
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 08:58:33 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

  Assert(it != method.end(), tk::ExceptType::FATAL,
         std::string("Cannot find parameter for MKLGaussianMethod \"") + m + "\"");

  return it->second;
}
