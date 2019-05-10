// *****************************************************************************
/*!
  \file      src/DiffEq/Velocity/ConfigureVelocity.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register velocity SDE into DiffEq factory
void registerVelocity( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the velocity SDE
std::vector< std::pair< std::string, std::string > >
infoVelocity( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureVelocity_h
