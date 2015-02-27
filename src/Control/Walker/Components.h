//******************************************************************************
/*!
  \file      src/Control/Walker/Components.h
  \author    J. Bakosi
  \date      Fri 27 Feb 2015 06:46:12 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
//******************************************************************************
#ifndef WalkerComponents_h
#define WalkerComponents_h

#include <Components.h>

namespace walker {
namespace ctr {

//! Number of components of systems of equations
using ncomps = tk::ctr::ncomponents<
  tag::dirichlet,    std::vector< tk::ctr::ncomp_type >,
  tag::gendir,       std::vector< tk::ctr::ncomp_type >,
  tag::wrightfisher, std::vector< tk::ctr::ncomp_type >,
  tag::diagou,       std::vector< tk::ctr::ncomp_type >,
  tag::ou,           std::vector< tk::ctr::ncomp_type >,
  tag::skewnormal,   std::vector< tk::ctr::ncomp_type >,
  tag::gamma,        std::vector< tk::ctr::ncomp_type >,
  tag::beta,         std::vector< tk::ctr::ncomp_type >,
  tag::nfracbeta,    std::vector< tk::ctr::ncomp_type >
>;

} // ctr::
} // walker::

#endif // WalkerComponents_h
