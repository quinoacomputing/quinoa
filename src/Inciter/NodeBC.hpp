// *****************************************************************************
/*!
  \file      src/Inciter/NodeBC.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions for nodal discretizations
  \details   Boundary conditions for nodal discretizations, such as continuous
    Galerkin finite elements.
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

} // inciter::

#endif // NodeBC_h
