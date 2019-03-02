// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for compressible flow PDE
  \details   Register and compile configuration for compressible flow PDE.
*/
// *****************************************************************************
#ifndef ConfigureCompFlow_h
#define ConfigureCompFlow_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.h"
#include "SystemComponents.h"
#include "Inciter/Options/PDE.h"

namespace inciter {

//! Register compressible flow PDEs into PDE factory
void
registerCompFlow( CGFactory& cf,
                  DGFactory& df,
                  std::set< ctr::PDEType >& cgt,
                  std::set< ctr::PDEType >& dgt );

//! Return information on the compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoCompFlow( std::map< ctr::PDEType, tk::ctr::ncomp_t >& cnt );

} // inciter::

#endif // ConfigureCompFlow_h
