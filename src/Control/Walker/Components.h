// *****************************************************************************
/*!
  \file      src/Control/Walker/Components.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during the control file parsing.
*/
// *****************************************************************************
#ifndef WalkerComponents_h
#define WalkerComponents_h

#include "SystemComponents.h"

namespace walker {
namespace ctr {

//! Number of components of systems of equations
using ncomps = tk::ctr::ncomponents<
   tag::dirichlet,       std::vector< tk::ctr::ncomp_type >
 , tag::mixdirichlet,    std::vector< tk::ctr::ncomp_type >
 , tag::gendir,          std::vector< tk::ctr::ncomp_type >
 , tag::wrightfisher,    std::vector< tk::ctr::ncomp_type >
 , tag::diagou,          std::vector< tk::ctr::ncomp_type >
 , tag::ou,              std::vector< tk::ctr::ncomp_type >
 , tag::skewnormal,      std::vector< tk::ctr::ncomp_type >
 , tag::gamma,           std::vector< tk::ctr::ncomp_type >
 , tag::beta,            std::vector< tk::ctr::ncomp_type >
 , tag::numfracbeta,     std::vector< tk::ctr::ncomp_type >
 , tag::massfracbeta,    std::vector< tk::ctr::ncomp_type >
 , tag::mixnumfracbeta,  std::vector< tk::ctr::ncomp_type >
 , tag::mixmassfracbeta, std::vector< tk::ctr::ncomp_type >
 , tag::velocity,        std::vector< tk::ctr::ncomp_type >
 , tag::position,        std::vector< tk::ctr::ncomp_type >
 , tag::dissipation,     std::vector< tk::ctr::ncomp_type >
>;

} // ctr::
} // walker::

#endif // WalkerComponents_h
