//******************************************************************************
/*!
  \file      src/Control/GammaParameters.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 07:15:00 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gamma mix model parameters
  \details   Gamma mix model parameters
*/
//******************************************************************************
#ifndef GammaParameters_h
#define GammaParameters_h

#include <vector>

namespace quinoa {

struct GammaParameters {
  real c1;
  real c2;
  real c3;
  real c4;

  //! Default constructor with default parameter values
  explicit GammaParameters() : c1(0.5), c2(0.73), c3(5.0), c4(0.25) {}
};

} // namespace quinoa

#endif // GammaParameters_h
