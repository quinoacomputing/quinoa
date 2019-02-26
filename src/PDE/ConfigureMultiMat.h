// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiMat.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Register and compile configuration for multi-material compressible
     flow PDE
  \details   Register and compile configuration for compressible multi-material
     flow PDE.
*/
// *****************************************************************************
#ifndef ConfigureMultiMat_h
#define ConfigureMultiMat_h

#include <set>
#include <map>
#include <vector>

#include "PDEFactory.h"
#include "SystemComponents.h"
#include "Inciter/Options/PDE.h"

namespace inciter {

//! Register compressible flow PDEs into PDE factory
void
registerMultiMat( DGFactory& df, std::set< ctr::PDEType >& dgt );

//! Return information on the multi-material compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoMultiMat( std::map< ctr::PDEType, tk::ctr::ncomp_type >& cnt );

} // inciter::

#endif // ConfigureMultiMat_h
