// *****************************************************************************
/*!
  \file      src/Mesh/CommMap.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Type definitions for communication maps
  \details   Type definitions for communication maps.
*/
// *****************************************************************************
#ifndef CommMap_h
#define CommMap_h

#include <unordered_set>
#include <unordered_map>
#include <map>

#include "Tags.hpp"
#include "UnsMesh.hpp"
#include "TaggedTuple.hpp"

namespace tk {

//! Node list type used for node communication map
using NodeSet = std::unordered_set< std::size_t >;

//! Edge list type used for edge communication map
using EdgeSet = UnsMesh::EdgeSet;

//! Node communication map type
//! \details Global mesh node IDs bordering the mesh chunk held by fellow
//!    chares associated to chare IDs
using NodeCommMap = std::unordered_map< int, NodeSet >;

//! Edge communication map type
//! \details Edge-end points with global mesh node IDs bordering the mesh chunk
//!   held by fellow chares associated to chare IDs
using EdgeCommMap = std::unordered_map< int, EdgeSet >;

//! Type list of all communication maps bundled as a tagged tuple
using AllCommMaps =
  TaggedTuple< brigand::list<
     tag::node,    NodeSet
   , tag::edge,    EdgeSet
  > >;

//! Communication map bundle type
//! \details All types of communication maps bundled and associated to chare IDs
using CommMaps = std::map< int, AllCommMaps >;

//! Decide if a node is not counted by a chare
bool slave( const NodeCommMap& map, std::size_t node, int chare );

//! Count the number of contributions to a node
tk::real count( const NodeCommMap& map, std::size_t node );

} // tk::

#endif // CommMap_h
