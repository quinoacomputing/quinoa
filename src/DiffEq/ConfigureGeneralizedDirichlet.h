// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureGeneralizedDirichlet.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the generalized Dirichlet SDE
  \details   Register and compile configuration on the generalized Dirichlet SDE.
*/
// *****************************************************************************
#ifndef ConfigureGeneralizedDirichlet_h
#define ConfigureGeneralizedDirichlet_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register generalized Dirichlet SDE into DiffEq factory
void registerGenDir( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the generlized Dirichlet SDE
std::vector< std::pair< std::string, std::string > >
infoGenDir( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureGenealizedDirichlet_h
