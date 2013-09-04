//******************************************************************************
/*!
  \file      src/Control/BetaParameters.h
  \author    J. Bakosi
  \date      Tue Sep  3 07:39:13 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model parameters
  \details   Beta mass model parameters
*/
//******************************************************************************
#ifndef BetaParameters_h
#define BetaParameters_h

namespace quinoa {

struct BetaParameters {
  real atwood;

  //! Default constructor with default parameter values
  explicit BetaParameters() : atwood(0.0) {}
};

} // namespace quinoa

#endif // BetaParameters_h
