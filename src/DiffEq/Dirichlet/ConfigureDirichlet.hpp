// *****************************************************************************
/*!
  \file      src/DiffEq/Dirichlet/ConfigureDirichlet.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the Dirichlet SDE
  \details   Register and compile configuration on the Dirichlet SDE.
*/
// *****************************************************************************
#ifndef ConfigureDirichlet_h
#define ConfigureDirichlet_h

#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register Dirichlet SDE into DiffEq factory
void registerDirichlet( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the Dirichlet SDE
std::vector< std::pair< std::string, std::string > >
infoDirichlet( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureDirichlet_h
