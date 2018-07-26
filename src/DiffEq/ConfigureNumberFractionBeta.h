// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureNumberFractionBeta.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the number fraction beta SDE
  \details   Register and compile configuration on the number fraction beta SDE.
*/
// *****************************************************************************
#ifndef ConfigureNumberFractionBeta_h
#define ConfigureNumberFractionBeta_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register beta SDE into DiffEq factory
void registerNumberFractionBeta( DiffEqFactory& f,
                                 std::set< ctr::DiffEqType >& t );

//! Return information on the beta SDE
std::vector< std::pair< std::string, std::string > >
infoNumberFractionBeta( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureNumberFractionBeta_h
