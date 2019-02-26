// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureMixNumberFractionBeta.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration on the mix number fraction beta
             SDE
  \details   Register and compile configuration on the mix number fraction beta
             SDE.
*/
// *****************************************************************************
#ifndef ConfigureMixNumberFractionBeta_h
#define ConfigureMixNumberFractionBeta_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register beta SDE into DiffEq factory
void registerMixNumberFractionBeta( DiffEqFactory& f,
                                    std::set< ctr::DiffEqType >& t );

//! Return information on the beta SDE
std::vector< std::pair< std::string, std::string > >
infoMixNumberFractionBeta( std::map< ctr::DiffEqType,
                                     tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureMixNumberFractionBeta_h
