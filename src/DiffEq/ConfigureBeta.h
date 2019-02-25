// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureBeta.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the beta SDE
  \details   Register and compile configuration on the beta SDE.
*/
// *****************************************************************************
#ifndef ConfigureBeta_h
#define ConfigureBeta_h

#include <set>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register beta SDE into DiffEq factory
void registerBeta( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the beta SDE
std::vector< std::pair< std::string, std::string > >
infoBeta( std::map< ctr::DiffEqType, tk::ctr::ncomp_t >& cnt );

} // walker::

#endif // ConfigureBeta_h
