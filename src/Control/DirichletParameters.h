//******************************************************************************
/*!
  \file      src/Control/DirichletParameters.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 07:13:23 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model parameters
  \details   Dirichlet mix model parameters
*/
//******************************************************************************
#ifndef DirichletParameters_h
#define DirichletParameters_h

#include <vector>

namespace quinoa {

struct DirichletParameters {
  std::vector<real> b;
  std::vector<real> S;
  std::vector<real> kappa;

  //! Default constructor with default parameter values
  explicit DirichletParameters() : b(), S(), kappa() {}
};

} // namespace quinoa

#endif // DirichletParameters_h
