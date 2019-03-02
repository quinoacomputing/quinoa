// *****************************************************************************
/*!
  \file      src/Control/Walker/Components.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Storage for number of components
  \details   Storage for number of components. This is part of the input deck
     stack and is thus populated during control file parsing.
*/
// *****************************************************************************
#ifndef WalkerComponents_h
#define WalkerComponents_h

#include "SystemComponents.h"

namespace walker {
namespace ctr {

//! Number of components storage for all systems of equations supported
using ncomps = tk::ctr::ncomponents<
   tag::dirichlet
 , tag::mixdirichlet
 , tag::gendir
 , tag::wrightfisher
 , tag::diagou
 , tag::ou
 , tag::skewnormal
 , tag::gamma
 , tag::beta
 , tag::numfracbeta
 , tag::massfracbeta
 , tag::mixnumfracbeta
 , tag::mixmassfracbeta
 , tag::velocity
 , tag::position
 , tag::dissipation
>;

} // ctr::
} // walker::

#endif // WalkerComponents_h
