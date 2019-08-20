// *****************************************************************************
/*!
  \file      src/DiffEq/Gamma/ConfigureGamma.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the gamma SDE
  \details   Register and compile configuration on the gamma SDE.
*/
// *****************************************************************************
#ifndef ConfigureGamma_h
#define ConfigureGamma_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register gamma SDE into DiffEq factory
void registerGamma( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the gamma SDE
std::vector< std::pair< std::string, std::string > >
infoGamma( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureGamma_h
