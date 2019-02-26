// *****************************************************************************
/*!
  \file      src/DiffEq/ConfigureSkewNormal.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration on the skew-normal SDE
  \details   Register and compile configuration on the skew-normal SDE.
*/
// *****************************************************************************
#ifndef ConfigureSkewNormal_h
#define ConfigureSkewNormal_h

#include <set>
#include <map>
#include <string>
#include <utility>

#include "DiffEqFactory.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Register skew-normal SDE into DiffEq factory
void registerSkewNormal( DiffEqFactory& f, std::set< ctr::DiffEqType >& t );

//! Return information on the skew-normal SDE
std::vector< std::pair< std::string, std::string > >
infoSkewNormal( std::map< ctr::DiffEqType, tk::ctr::ncomp_type >& cnt );

} // walker::

#endif // ConfigureSkewNormal_h
