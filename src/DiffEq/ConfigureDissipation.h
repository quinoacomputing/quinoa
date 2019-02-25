// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureDissipation.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the dissipation SDE
  \details   Register and compile configuration on the dissipation SDE.
*/
// *****************************************************************************
#ifndef ConfigureDissipation_h
#define ConfigureDissipation_h

#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register dissipation SDE into DiffEq factory
void registerDissipation( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the dissipation SDE
std::vector< std::pair< std::string, std::string > >
infoDissipation( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureDissipation_h
