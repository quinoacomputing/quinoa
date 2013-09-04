//******************************************************************************
/*!
  \file      src/Control/GLMParameters.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 07:15:35 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydro model parameters
  \details   Generalized Langevin hydro model parameters
*/
//******************************************************************************
#ifndef GLMParameters_h
#define GLMParameters_h

namespace quinoa {

struct GLMParameters {
  real c0;

  //! Default constructor with default parameter values
  explicit GLMParameters() : c0(2.1) {}
};

} // namespace quinoa

#endif // GLMParameters_h
