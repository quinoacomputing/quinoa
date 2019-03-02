// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureDiagOrnsteinUhlenbeck.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the diagonal
             Ornstein-Uhlenbeck SDE
  \details   Register and compile configuration on the diagonal
             Ornstein-Uhlenbeck SDE.
*/
// *****************************************************************************
#ifndef ConfigureDiagOrnsteinUhlenbeck_h
#define ConfigureDiagOrnsteinUhlenbeck_h

#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register diagonal Ornstein-Uhlenbeck SDE into DiffEq factory
void registerDiagOrnsteinUhlenbeck( DiffEqFactory& f,
                                    std::set< ctr::DiffEqType >& t );

//! Return information on the diagonal Ornstein-Uhlenbeck SDE
std::vector< std::pair< std::string, std::string > >
infoDiagOrnsteinUhlenbeck( std::map< ctr::DiffEqType,
                                     tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureDiagOrnsteinUhlenbeck_h
