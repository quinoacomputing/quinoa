// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureMixDirichlet.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Register and compile configuration on the MixDirichlet SDE
  \details   Register and compile configuration on the MixDirichlet SDE.
*/
// *****************************************************************************
#ifndef ConfigureMixDirichlet_h
#define ConfigureMixDirichlet_h

#include <set>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register MixDirichlet SDE into DiffEq factory
void registerMixDirichlet( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the MixDirichlet SDE
std::vector< std::pair< std::string, std::string > >
infoMixDirichlet( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureMixDirichlet_h
