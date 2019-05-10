// *****************************************************************************
/*!
  \file      src/DiffEq/OrnsteinUhlenbeck/ConfigureOrnsteinUhlenbeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register Ornstein-Uhlenbeck SDE into DiffEq factory
void
registerOrnsteinUhlenbeck( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the Ornstein-Uhlenbeck SDE
std::vector< std::pair< std::string, std::string > >
infoOrnsteinUhlenbeck( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureOrnsteinUhlenbeck_h
