//******************************************************************************
/*!
  \file      src/Control/Walker/Components.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 07:24:57 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components
*/
//******************************************************************************
#ifndef WalkerComponents_h
#define WalkerComponents_h

#include <Components.h>

namespace walker {
namespace ctr {

//! Number of components of models and equations
using ncomps = tk::ctr::ncomponents<
  tag::dirichlet,    std::vector< unsigned int >, //!< Dirichlet SDEs
  tag::gendir,       std::vector< unsigned int >, //!< Generalized Dirichlet
  tag::wrightfisher, std::vector< unsigned int >, //!< Wright-Fisher SDEs
  tag::diagou,       std::vector< unsigned int >, //!< Diagonal OU SDEs
  tag::ou,           std::vector< unsigned int >, //!< Ornstein-Uhlenbeck SDEs
  tag::skewnormal,   std::vector< unsigned int >, //!< Skew-normal SDEs
  tag::gamma,        std::vector< unsigned int >, //!< Gamma SDEs
  tag::beta,         std::vector< unsigned int >  //!< Beta SDEs
>;

} // ctr::
} // walker::

#endif // WalkerComponents_h
