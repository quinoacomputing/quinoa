// *****************************************************************************
/*!
  \file      src/DiffEq/Beta/ConfigureMixMassFractionBeta.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the mix mass fraction beta
             SDE
  \details   Register and compile configuration on the mix mass fraction beta
             SDE.
*/
// *****************************************************************************
#ifndef ConfigureMixMassFractionBeta_h
#define ConfigureMixMassFractionBeta_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.hpp"
#include "Walker/Options/DiffEq.hpp"

namespace walker {

//! Register beta SDE into DiffEq factory
void registerMixMassFractionBeta( DiffEqFactory& f,
                                  std::set< ctr::DiffEqType >& t );

//! Return information on the beta SDE
std::vector< std::pair< std::string, std::string > >
infoMixMassFractionBeta( std::map< ctr::DiffEqType,
                                     tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureMixMassFractionBeta_h
