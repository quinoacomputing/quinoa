// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureOrnsteinUhlenbeck.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the Ornstein-Uhlenbeck SDE
  \details   Register and compile configuration on the Ornstein-Uhlenbeck SDE.
*/
// *****************************************************************************
#ifndef ConfigureOrnsteinUhlenbeck_h
#define ConfigureOrnsteinUhlenbeck_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register Ornstein-Uhlenbeck SDE into DiffEq factory
void
registerOrnsteinUhlenbeck( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the Ornstein-Uhlenbeck SDE
std::vector< std::pair< std::string, std::string > >
infoOrnsteinUhlenbeck( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureOrnsteinUhlenbeck_h
