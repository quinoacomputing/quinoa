// *****************************************************************************
/*!
  \file      src/Inciter/NodeBC.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Boundary conditions for nodal discretizations
  \details   Boundary conditions for nodal discretizations, such as continuous
    Galerkin finite elements.
*/
// *****************************************************************************

#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "NodeBC.hpp"
#include "CGPDE.hpp"
#include "Fields.hpp"
#include "Vector.hpp"

namespace inciter {

extern std::vector< CGPDE > g_cgpde;

std::unordered_map< std::size_t, std::vector< std::pair< bool, tk::real > > >
match( std::size_t meshid,
       [[maybe_unused]] tk::ncomp_t ncomp,
       tk::real t,
       tk::real dt,
       const std::vector< tk::real >& tp,
       const std::vector< tk::real >& dtp,
       const tk::UnsMesh::Coords& coord,
       const std::unordered_map< std::size_t, std::size_t >& lid,
       const std::map< int, std::vector< std::size_t > >& bnode,
       bool increment )
// *****************************************************************************
//  Match user-specified boundary conditions at nodes for side sets
//! \param[in] ncomp Number of scalar components in PDE system
//! \param[in] t Physical time at which to query boundary conditions
//! \param[in] dt Time step size (for querying BC increments in time)
//! \param[in] tp Physical time for each mesh node
//! \param[in] dtp Time step size for each mesh node
//! \param[in] coord Mesh node coordinates
//! \param[in] lid Local node IDs associated to local node IDs
//! \param[in] bnode Map storing global mesh node IDs mapped to side set ids
//! \param[in] increment If true, evaluate the solution increment between
//!   t and t+dt for Dirichlet BCs. If false, evlauate the solution instead.
//! \return Vector of pairs of bool and boundary condition value associated to
//!   local mesh node IDs at which the user has set Dirichlet boundary
//!   conditions for all systems of PDEs integrated. The bool indicates whether
//!   the BC is set at the node for that component: if true, the real value is
//!   the increment (from t to dt) in (or the value of) the BC specified for a
//!   component.
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets (i.e., a
//!   list of elements). This function queries Dirichlet boundary condition
//!   values from all PDEs in the multiple systems of PDEs integrated at the
//!   node lists associated to side set IDs, given by bnode. Each
//!   PDE system returns a BC data structure. Note that the BC mesh nodes that
//!   this function results in (stored in dirbc) only contains those nodes that
//!   are supplied via bnode. i.e., in parallel only a part of the mesh is
//!   worked on.
// *****************************************************************************
{
  using inciter::g_cgpde;

  // Vector of pairs of bool and boundary condition value associated to mesh
  // node IDs at which the user has set Dirichlet BCs for all PDEs integrated.
  std::unordered_map< std::size_t,
    std::vector< std::pair< bool, tk::real > > > dirbc;

  // Details for the algorithm below: PDE::dirbc() returns a new map that
  // associates a vector of pairs associated to local node IDs. (The pair is a
  // pair of bool and real value, the former is the fact that the BC is to be
  // set while the latter is the value if it is to be set). The length of this
  // NodeBC vector, returning from each system of PDEs equals to the number of
  // scalar components the given PDE integrates. Here we contatenate this map
  // for all PDEs being integrated. If there are multiple BCs set at a mesh node
  // (dirbc::key), either because (1) in the same PDE system the user prescribed
  // BCs on side sets that share nodes or (2) because more than a single PDE
  // system assigns BCs to a given node (on different variables), the NodeBC
  // vector must be correctly stored. "Correctly" here means that the size of
  // the NodeBC vectors must all be the same and equal to the sum of all scalar
  // components integrated by all PDE systems integrated. Example: single-phase
  // compressible flow (density, momentum, energy = 5) + transported scalars of
  // 10 variables -> NodeBC vector length = 15. Note that in case (1) above a
  // new node encountered must "overwrite" the already existing space for the
  // NodeBC vector. "Overwrite" here means that it should keep the existing BCs
  // and add the new ones yielding the union the two prescription for BCs but in
  // the same space that already exist in the NodeBC vector. In case (2),
  // however, the NodeBC pairs must go to the location in the vector assigned to
  // the given PDE system, i.e., using the above example BCs for the 10 (or
  // less) scalars should go in the positions starting at 5, leaving the first 5
  // false, indicating no BCs for the flow variables.
  //
  // Note that the logic described above is only partially implemented at this
  // point. What works is the correct insertion of multiple BCs for nodes shared
  // among multiple side sets, e.g., corners, originating from the same PDE
  // system. What is not yet implemented is the case when there are no BCs set
  // for flow variables but there are BCs for transport, the else branch below
  // will incorrectly NOT skip the space for the flow variables. In other words,
  // this only works for a single PDE system and a sytem of systems. This
  // machinery is only tested with a single system of PDEs at this point.
  //
  // When a particular node belongs to two or more side sets with different BCs,
  // there is an ambiguity as to which of the multiple BCs should be applied to
  // the node. This issue is described in case (1) above. In the current
  // implementation, every side set applies the BC to the common node in
  // question, successively overwriting the BC applied by the previous side set.
  // Effectively, the BC corresponding to the last side set ID is applied to the
  // common node. Since bnode is an ordered map, the side set with a larger
  // id wins if a node belongs to multiple side sets.

  // Lambda to convert global to local node ids of a list of nodes
  auto local = [ &lid ]( const std::vector< std::size_t >& gnodes ){
    std::vector< std::size_t > lnodes( gnodes.size() );
    for (std::size_t i=0; i<gnodes.size(); ++i)
      lnodes[i] = tk::cref_find( lid, gnodes[i] );
    return lnodes;
  };

  // Query Dirichlet BCs for all PDEs integrated and assign to nodes
  for (const auto& s : bnode) {     // for all side sets passed in
    std::size_t c = 0;
    auto l = local(s.second);   // generate local node ids on side set
    // query Dirichlet BCs at nodes of this side set
    auto eqbc =
      g_cgpde[meshid].dirbc( t, dt, tp, dtp, {s.first,l}, coord, increment );
    for (const auto& n : eqbc) {
      auto id = n.first;                      // BC node ID
      const auto& bcs = n.second;             // BCs
      auto& nodebc = dirbc[ id ];     // BCs to be set for node
      if (nodebc.size() > c) {        // node already has BCs from this PDE
        Assert( nodebc.size() == c+bcs.size(), "Size mismatch" );
        for (std::size_t i=0; i<bcs.size(); i++) {
          if (bcs[i].first) nodebc[c+i] = bcs[i];
        }
      } else {        // node does not yet have BCs from this PDE
        // This branch needs to be completed for system of systems of PDEs.
        // See note above.
        nodebc.insert( end(nodebc), begin(bcs), end(bcs) );
      }
    }
    if (!eqbc.empty()) c += eqbc.cbegin()->second.size();
  }

  // Verify the size of each NodeBC vectors. They must have the same lengths and
  // equal to the total number of scalar components for all systems of PDEs
  // integrated.
  Assert( std::all_of( begin(dirbc), end(dirbc),
            [ ncomp ]( const auto& n ){ return n.second.size() == ncomp; } ),
          "Size of NodeBC vector incorrect" );
 
  return dirbc;
}

} // inciter::
