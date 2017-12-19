// *****************************************************************************
/*!
  \file      src/Inciter/BoundaryConditions.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Data and functionality working on boundary conditions
  \details   Data and functionality working on boundary conditions.
*/
// *****************************************************************************

#include <vector>
#include <map>
#include <unordered_map>

#include "BoundaryConditions.h"

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

#include "NoWarning/boundaryconditions.def.h"
