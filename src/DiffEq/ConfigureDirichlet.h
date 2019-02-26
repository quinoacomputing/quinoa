// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureDirichlet.h
  \copyright 2016-2018, Triad National Security, LLC.
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

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register Dirichlet SDE into DiffEq factory
void registerDirichlet( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the Dirichlet SDE
std::vector< std::pair< std::string, std::string > >
infoDirichlet( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureDirichlet_h
