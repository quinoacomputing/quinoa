// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureWrightFisher.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration on the Wright-Fisher SDE
  \details   Register and compile configuration on the Wright-Fisher SDE.
*/
// *****************************************************************************
#ifndef ConfigureWrightFisher_h
#define ConfigureWrightFisher_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register Wright-Fisher SDE into DiffEq factory
void registerWrightFisher( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the Wright-Fisher SDE
std::vector< std::pair< std::string, std::string > >
infoWrightFisher( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureWrightFisher_h
