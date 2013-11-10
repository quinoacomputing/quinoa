//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MKLGaussianMethod.C
  \author    J. Bakosi
  \date      Sat 09 Nov 2013 06:08:23 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options
  \details   Quinoa's random number generator options
*/
//******************************************************************************

#include <Quinoa/Options/MKLGaussianMethod.h>

using namespace quinoa::ctr;

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
