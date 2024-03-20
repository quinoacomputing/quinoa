// *****************************************************************************
/*!
  \file      src/Inciter/NodeBC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions for nodal discretizations
  \details   Boundary conditions for nodal discretizations, such as continuous
    Galerkin finite elements, e.g., DiagCG.
*/
// *****************************************************************************
#ifndef NodeBC_h
#define NodeBC_h

#include <vector>
#include <map>
#include <unordered_map>

#include "Types.hpp"
#include "UnsMesh.hpp"
#include "Fields.hpp"

namespace inciter {

//! Match user-specified boundary conditions at nodes for side sets
std::unordered_map< std::size_t, std::vector< std::pair< bool, tk::real > > >
match( std::size_t meshid,
       tk::ncomp_t ncomp,
       tk::real t,
       tk::real dt,
       const std::vector< tk::real >& tp,
       const std::vector< tk::real >& dtp,
       const tk::UnsMesh::Coords& coord,
       const std::unordered_map< std::size_t, std::size_t >& lid,
       const std::map< int, std::vector< std::size_t > >& sidenodes,
       bool increment );

//! \brief Verify that the change in the solution at those nodes where
//!   Dirichlet boundary conditions are set is exactly the amount the BCs
//!   prescribe
bool
correctBC( const tk::Fields& a,
           const tk::Fields& dul,
           const std::unordered_map< std::size_t,
                   std::vector< std::pair< bool, tk::real > > >& bc );

} // inciter::

#endif // NodeBC_h
