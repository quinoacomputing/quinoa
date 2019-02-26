// *****************************************************************************
/*!
  \file      src/PDE/ConfigureCompFlow.h
  \copyright 2016-2018, Triad National Security, LLC.
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
infoCompFlow( std::map< ctr::PDEType, tk::ctr::ncomp_type >& cnt );

} // inciter::

#endif // ConfigureCompFlow_h
