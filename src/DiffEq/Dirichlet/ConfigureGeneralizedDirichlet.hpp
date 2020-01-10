// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/ConfigureGeneralizedDirichlet.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register generalized Dirichlet SDE into DiffEq factory
void registerGenDir( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the generlized Dirichlet SDE
std::vector< std::pair< std::string, std::string > >
infoGenDir( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureGenealizedDirichlet_h
