// *****************************************************************************
/*!
  \file      src/PDE/ConfigureTransport.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration on the Transport PDE
  \details   Register and compile configuration on the Transport PDE.
*/
// *****************************************************************************
#ifndef ConfigureTransport_h
#define ConfigureTransport_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.h"
#include "SystemComponents.h"
#include "Inciter/Options/PDE.h"

namespace inciter {

//! Register transport PDEs into PDE factory
void
registerTransport( CGFactory& cf,
                   DGFactory& df,
                   std::set< ctr::PDEType >& cgt,
                   std::set< ctr::PDEType >& dgt );

//! Return information on the transport PDE
std::vector< std::pair< std::string, std::string > >
infoTransport( std::map< ctr::PDEType, tk::ctr::ncomp_type >& cnt );

} // inciter::

#endif // ConfigureTransport_h
