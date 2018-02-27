// *****************************************************************************
/*!
  \file      src/Inciter/BoundaryConditions.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Data and functionality working on boundary conditions
  \details   Data and functionality working on boundary conditions.
*/
// *****************************************************************************

#include <vector>
#include <map>
#include <unordered_map>

#include "BoundaryConditions.h"
#include "SystemComponents.h"
#include "CGPDE.h"

namespace inciter {

extern std::vector< CGPDE > g_cgpde;

}

using inciter::BoundaryConditions;

BoundaryConditions::BoundaryConditions(
  const std::map< int, std::vector< std::size_t > >& sidenodes )
  : m_sideFileNodes( sidenodes )
// *****************************************************************************
//  Constructor
//! \param[in] ss File mesh node IDs mapped to side set ids
// *****************************************************************************
{
}

std::map< int, std::vector< std::size_t > >
BoundaryConditions::sideNodes(
  const std::unordered_map< std::size_t, std::size_t >& filenodes,
  const std::unordered_map< std::size_t, std::size_t >& lid )
// *****************************************************************************
// Create map that assigns the local mesh node IDs mapped to side set ids
//! \param[in] filenodes Map associating file node IDs to local node IDs
//! \param[in] lid Local node IDs associated to global node IDs
//! \return Map that assigns the local mesh node IDs mapped to side set ids,
//!   storing only those nodes for a given side set that are part of our chunk
//!   of the mesh (based on a search in filenodes)
// *****************************************************************************
{
  // First generate map associating local node IDs to file node IDs. We invert
  // the map that associates file node IDs to local node IDs for the purpose of
  // enabling efficient searches of file node IDs.
  std::unordered_map< std::size_t, std::size_t > localnodes;
  for (const auto& i : filenodes) {
    auto n = tk::cref_find( lid, i.first );
    Assert( n < lid.size(),
            "Local IDs must be lower than the local number of grid points" );
    localnodes[ i.second ] = n;
  }

  // Create map that assigns the local mesh node IDs mapped to side set ids
  std::map< int, std::vector< std::size_t > > sidenodes;
  for (const auto& s : m_sideFileNodes) {
    auto& n = sidenodes[ s.first ];
    for (auto o : s.second) {
      auto it = localnodes.find( o );
      if (it != end(localnodes))
        n.push_back( it->second );
    }
  }

  return sidenodes;
}

std::unordered_map< std::size_t,
  std::vector< std::pair< bool, tk::real > > >
BoundaryConditions::match( tk::ctr::ncomp_type ncomp,
                           tk::real t,
                           tk::real dt,
                           const tk::UnsMesh::Coords& coord,
                           const std::vector< std::size_t > gid,
                           const std::map< int, std::vector< std::size_t > >&
                             sidenodes )
// *****************************************************************************
//  Query and matchto user-specified boundary conditions to side sets
//! \param[in] ncomp Number of scalar components in PDE system
//! \param[in] t Physical time at which to query boundary conditions
//! \param[in] dt Time step size (for querying BC increments in time)
//! \param[in] coord Mesh node coordinates
//! \param[in] gid Global node IDs
//! \param[in] sidenodes Map storing local mesh node IDs mapped to side set ids
//! \return Vector of pairs of bool and boundary condition value associated to mesh
//!   node IDs at which the user has set Dirichlet boundary conditions for all
//!   PDEs integrated. The bool indicates whether the BC is set at the node for
//!   that component the if true, the real value is the increment (from t to dt)
//!   in the BC specified for a component.
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets (i.e., a
//!   list of elements). This function queries Dirichlet boundary condition
//!   values from all PDEs in the system of PDEs integrated at the node lists
//!   associated to side set IDs (previously read from file). As a response to
//!   this query, each PDE system returns a BC data structure which is then
//!   sent to the linear system solver which needs to know about this to apply
//!   BCs before a linear solve. Note that the BC mesh nodes that this function
//!   results in, stored in dirbc and sent to the linear system solver, only
//!   contains those nodes that this chare contributes to, i.e., it does not
//!   contain those BC nodes at which other chares enforce Dirichlet BCs. The
//!   linear system solver then collects these and communicates to other PEs so
//!   that BC data held in Solver::m_bc are the same on all PEs.
// *****************************************************************************
{
  // Vector of pairs of bool and boundary condition value associated to mesh
  // node IDs at which the user has set Dirichlet boundary conditions for all
  // PDEs integrated. NodeBC = value_type.
  std::unordered_map< std::size_t,
    std::vector< std::pair< bool, tk::real > > > dirbc;

  // Query Dirichlet boundary conditions for all PDEs integrated and assign to
  // nodes. This is where the individual system of PDEs are queried for boundary
  // conditions. The outer loop goes through all sides sets that exists in the
  // input file and passes the map's value_type (a pair of the side set id and a
  // vector of local node IDs) to PDE::dirbc(). PDE::dirbc() returns a new map
  // that associates a vector of pairs associated to local node IDs. (The pair
  // is a pair of bool and real value, the former is the fact that the BC is to
  // be set while the latter is the value if it is to be set). The length of
  // this NodeBC vector, returning from each system of PDEs equals to the number
  // of scalar components the given PDE integrates. Here then we contatenate
  // this map for all PDEs integrated. If there are multiple BCs set at a mesh
  // node (dirbc::key), either because (1) in the same PDE system the user
  // prescribed BCs on side sets that share nodes or (2) because more than a
  // single PDE system assigns BCs to a given node (on different variables), the
  // NodeBC vector must be correctly stored. "Correctly" here means that the
  // size of the NodeBC vectors must all be the same and qual to the sum of all
  // scalar components integrated by all PDE systems integrated. Example:
  // single-phase compressible flow (density, momentum, energy = 5) +
  // transported scalars of 10 variables -> NodeBC vector length = 15. Note that
  // in case (1) above a new node encountered must "overwrite" the already
  // existing space for the NodeBC vector. "Overwrite" here means that it should
  // keep the existing BCs and add the new ones yielding the union the two
  // prescription for BCs but in the same space that already exist in the NodeBC
  // vector. In case (2), however, the NodeBC pairs must go to the location in
  // the vector assigned to the given PDE system, i.e., using the above example
  // BCs for the 10 (or less) scalars should go in the positions starting at 5,
  // leaving the first 5 false, indicating no BCs for the flow variables.
  //
  // TODO: Note that the logic described above is only partially implemented at
  // this point. What works is the correct insertion of multiple BCs for nodes
  // shared among multiple side sets, e.g., corners, originating from the same
  // PDE system. What is not yet implemented is the case when there are no BCs
  // set for flow variables but there are BCs for transport, the else branch
  // below will incorrectly NOT skip the space for the flow variables. In other
  // words, this only works for a single PDE system and a sytem of systems. This
  // machinery is only tested with a single system of PDEs at this point.

  for (const auto& s : sidenodes) {
    std::size_t c = 0;
    for (std::size_t eq=0; eq<g_cgpde.size(); ++eq) {
      auto eqbc = g_cgpde[eq].dirbc( t, dt, s, coord );
      for (const auto& n : eqbc) {
        auto id = n.first;                      // BC node ID
        const auto& bcs = n.second;             // BCs
        auto& nodebc = dirbc[ gid[id] ];        // BCs to be set for node
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
  }

  // Verify the size of each NodeBC vectors. They must have the same lengths and
  // equal to the total number of scalar components for all systems of PDEs
  // integrated. This is intentional, because this way the linear system solver
  // does not have to (and does not) know about individual equation systems.
  // This entire loop is optimized away in RELEASE mode, thus the IGNOREs to
  // silence compiler warnings.
  for (const auto& n : dirbc) {
    IGNORE(n);
    Assert( n.second.size() == ncomp, "Size of NodeBC vector incorrect" );
  }
  IGNORE(ncomp);

  return dirbc;
}

#include "NoWarning/boundaryconditions.def.h"
