// *****************************************************************************
/*!
  \file      src/Inciter/NodeBC.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Boundary conditions for nodal discretizations
  \details   Boundary conditions for nodal discretizations, such as continuous
    Galerkin finite elements, e.g., MatCG, DiagCG.
*/
// *****************************************************************************
#ifndef NodeBC_h
#define NodeBC_h

#include <vector>
#include <map>
#include <unordered_map>

#include "SystemComponents.h"
#include "UnsMesh.h"
#include "Fields.h"

namespace inciter {

//! Match user-specified boundary conditions at nodes for side sets
std::unordered_map< std::size_t, std::vector< std::pair< bool, tk::real > > >
match( tk::ctr::ncomp_type ncomp,
       tk::real t,
       tk::real dt,
       const tk::UnsMesh::Coords& coord,
       const std::vector< std::size_t > gid,
       const std::unordered_map< std::size_t, std::size_t >& lid,
       const std::map< int, std::vector< std::size_t > >& sidenodes );

//! \brief Verify that the change in the solution at those nodes where
//!   Dirichlet boundary conditions are set is exactly the amount the BCs
//!   prescribe
bool
correctBC( const tk::Fields& a,
           const tk::Fields& dul,
           const std::map< int, std::vector< std::size_t > >& bnode,
           const std::unordered_map< std::size_t,
                   std::vector< std::pair< bool, tk::real > > >& bc,
           const std::unordered_map< std::size_t, std::size_t >& lid );

} // inciter::

#endif // NodeBC_h
