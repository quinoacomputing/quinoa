// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureVelocity.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the velocity SDE
  \details   Register and compile configuration on the velocity SDE.
*/
// *****************************************************************************
#ifndef ConfigureVelocity_h
#define ConfigureVelocity_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register velocity SDE into DiffEq factory
void registerVelocity( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the velocity SDE
std::vector< std::pair< std::string, std::string > >
infoVelocity( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureVelocity_h
