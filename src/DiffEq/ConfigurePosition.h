// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigurePosition.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the position SDE
  \details   Register and compile configuration on the position SDE.
*/
// *****************************************************************************
#ifndef ConfigurePosition_h
#define ConfigurePosition_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register position SDE into DiffEq factory
void registerPosition( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the position SDE
std::vector< std::pair< std::string, std::string > >
infoPosition( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigurePosition_h
