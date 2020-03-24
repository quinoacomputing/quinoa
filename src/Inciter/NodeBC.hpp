// *****************************************************************************
/*!
  \file      src/Inciter/NodeBC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

#include "SystemComponents.hpp"
#include "UnsMesh.hpp"
#include "Fields.hpp"

namespace inciter {

//! Match user-specified boundary conditions at nodes for side sets
std::unordered_map< std::size_t, std::vector< std::pair< bool, tk::real > > >
match( tk::ctr::ncomp_t ncomp,
       tk::real t,
       tk::real dt,
       const tk::UnsMesh::Coords& coord,
       const std::unordered_map< std::size_t, std::size_t >& lid,
       const std::map< int, std::vector< std::size_t > >& sidenodes );

//! \brief Verify that the change in the solution at those nodes where
//!   Dirichlet boundary conditions are set is exactly the amount the BCs
//!   prescribe
bool
correctBC( const tk::Fields& a,
           const tk::Fields& dul,
           const std::unordered_map< std::size_t,
                   std::vector< std::pair< bool, tk::real > > >& bc );

//! Decide if node is a stagnation point
bool
stagPoint( const std::array< tk::real, 3 >& p,
           const std::tuple< std::vector< tk::real >,
                             std::vector< tk::real > >& stag );

} // inciter::

#endif // NodeBC_h
