//******************************************************************************
/*!
  \file      src/Control/GenDirichletParameters.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 07:14:29 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Dirichlet mix model parameters
  \details   Generalized Dirichlet mix model parameters
*/
//******************************************************************************
#ifndef GenDirichletParameters_h
#define GenDirichletParameters_h

#include <vector>

namespace quinoa {

struct GenDirichletParameters {
  std::vector<real> b;
  std::vector<real> S;
  std::vector<real> kappa;
  std::vector<real> c;

  //! Default constructor with default parameter values
  explicit GenDirichletParameters() : b(), S(), kappa(), c() {}
};

} // namespace quinoa

#endif // GenDirichletParameters_h
