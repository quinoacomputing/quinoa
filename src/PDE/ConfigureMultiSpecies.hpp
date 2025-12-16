// *****************************************************************************
/*!
  \file      src/PDE/ConfigureMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Register and compile configuration for multi-species compressible
     flow PDE
*/
// *****************************************************************************
#ifndef ConfigureMultiSpecies_h
#define ConfigureMultiSpecies_h

#include "PDEFactory.hpp"
#include "Inciter/Options/PDE.hpp"

namespace inciter {

//! Register compressible flow PDEs into PDE factory
void
registerMultiSpecies( DGFactory& df, FVFactory& ff,
  std::set< ctr::PDEType >& fvt, std::set< ctr::PDEType >& dgt );

//! Return information on the multi-species compressible flow PDE
std::vector< std::pair< std::string, std::string > >
infoMultiSpecies( std::map< ctr::PDEType, tk::ncomp_t >& cnt );

} // inciter::

#endif // ConfigureMultiSpecies_h
