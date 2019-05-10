// *****************************************************************************
/*!
  \file      src/PDE/ConfigureTransport.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration on the Transport PDE
  \details   Register and compile configuration on the Transport PDE.
*/
// *****************************************************************************
#ifndef ConfigureTransport_h
#define ConfigureTransport_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.hpp"
#include "SystemComponents.hpp"
#include "Inciter/Options/PDE.hpp"

namespace inciter {

//! Register transport PDEs into PDE factory
void
registerTransport( CGFactory& cf,
                   DGFactory& df,
                   std::set< ctr::PDEType >& cgt,
                   std::set< ctr::PDEType >& dgt );

//! Return information on the transport PDE
std::vector< std::pair< std::string, std::string > >
infoTransport( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt );

} // inciter::

#endif // ConfigureTransport_h
