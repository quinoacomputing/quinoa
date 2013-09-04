//******************************************************************************
/*!
  \file      src/Control/SLMParameters.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 07:15:06 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydro model parameters
  \details   Simplified Langevin hydro model parameters
*/
//******************************************************************************
#ifndef SLMParameters_h
#define SLMParameters_h

namespace quinoa {

struct SLMParameters {
  real c0;

  //! Default constructor with default parameter values
  explicit SLMParameters() : c0(2.1) {}
};

} // namespace quinoa

#endif // SLMParameters_h
