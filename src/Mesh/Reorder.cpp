// *****************************************************************************
/*!
  \file      src/Mesh/Reorder.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh reordering routines for unstructured meshes
  \details   Mesh reordering routines for unstructured meshes.
*/
// *****************************************************************************

#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <map>
#include <tuple>
#include <cstddef>

#include "Reorder.hpp"
#include "Exception.hpp"
#include "ContainerUtil.hpp"
#include "Vector.hpp"

namespace tk {

std::size_t
shiftToZero( std::vector< std::size_t >& inpoel )
// *****************************************************************************
//  Shift node IDs to start with zero in element connectivity
//! \param[inout] inpoel Inteconnectivity of points and elements
//! \return Amount shifted
//! \details This function implements a simple reordering of the node ids of the
//!   element connectivity in inpoel by shifting the node ids so that the
//!   smallest is zero.
//! \note It is okay to call this function with an empty container; it will
//!    simply return without throwing an exception.
// *****************************************************************************
{
  if (inpoel.empty()) return 0;

  // find smallest node id
  auto minId = *std::min_element( begin(inpoel), end(inpoel) );

  // shift node ids to start from zero
  // cppcheck-suppress useStlAlgorithm
  for (auto& n : inpoel) n -= minId;

  return minId;
}

void
remap( std::vector< std::size_t >& ids, const std::vector< std::size_t >& map )
// *****************************************************************************
//  Apply new maping to vector of indices
//! \param[inout] ids Vector of integer IDs to remap
//! \param[in] map Array of indices creating a new order
//! \details This function applies a mapping (reordering) to the integer IDs
//!   passed in using the map passed in. The mapping is expressed between the
//!   array index and its value. The function overwrites every value, i, of
//!   vector ids with map[i].
//! \note The sizes of ids and map need not equal. Only the maximum index in ids
//!   must be lower than the size of map.
//! \note It is okay to call this function with either of the containers empty;
//!   it will simply return without throwing an exception.
// *****************************************************************************
{
  if (ids.empty() || map.empty()) return;

  Assert( *max_element( begin(ids), end(ids) ) < map.size(),
          "Indexing out of bounds" );

  // remap integer IDs in vector ids
  // cppcheck-suppress useStlAlgorithm
  for (auto& i : ids) i = map[i];
}

void
remap( std::vector< tk::real >& r, const std::vector< std::size_t >& map )
// *****************************************************************************
//  Apply new maping to vector of real numbers
//! \param[inout] r Vector of real numbers to remap
//! \param[in] map Array of indices creating a new order
//! \details This function applies a mapping (reordering) to the real values
//!   passed in using the map passed in. The mapping is expressed between the
//!   array index and its value. The function moves every value r[i] to
//!   r[ map[i] ].
//! \note The sizes of r and map must be equal and the maximum index in map must
//!   be lower than the size of map.
//! \note It is okay to call this function with either of the containers empty;
//!   it will simply return without throwing an exception.
// *****************************************************************************
{
  if (r.empty() || map.empty()) return;

  Assert( r.size() == map.size(), "Size mismatch" );
  Assert( *max_element( begin(map), end(map) ) < map.size(),
          "Indexing out of bounds" );

  // remap real numbers in vector
  auto m = r;
  for (std::size_t i=0; i<map.size(); ++i) r[ map[i] ] = m[ i ];
}

std::vector< std::size_t >
remap( const std::vector< std::size_t >& ids,
       const std::vector< std::size_t >& map )
// *****************************************************************************
//  Create remapped vector of indices using a vector
//! \param[in] ids Vector of integer IDs to remap
//! \param[in] map Array of indices creating a new order
//! \return Remapped vector of ids
//! \details This function applies a mapping (reordering) to the integer IDs
//!   passed in using the map passed in. The mapping is expressed between the
//!   array index and its value. The function creates and returns a new container
//!   with remapped ids of identical size of the origin ids container.
//! \note The sizes of ids and map must be equal and the maximum index in map
//!   must be lower than the size of map.
//! \note It is okay to call this function with either of the containers empty;
//!   if ids is empty, it returns an empty container; if map is empty, it will
//!   return the original container.
// *****************************************************************************
{
  if (ids.empty()) return {};
  if (map.empty()) return ids;

  Assert( *max_element( begin(ids), end(ids) ) < map.size(),
          "Indexing out of bounds" );

  // in terms of the in-place remap of a vector usinga vector
  auto newids = ids;
  remap( newids, map );

  return newids;
}

void
remap( std::vector< std::size_t >& ids,
       const std::unordered_map< std::size_t, std::size_t >& map )
// *****************************************************************************
//  In-place remap vector of indices using a map
//! \param[in] ids Vector of integer IDs to remap
//! \param[in] map Hash-map of key->value creating a new order
//! \details This function applies a mapping (reordering) to the integer IDs
//!   passed in using the map passed in. The mapping is expressed as a hash-map
//!   of key->value pairs, where the key is the original and the value is the
//!   new ids of the mapping. The function overwrites the ids container with the
//!   remapped ids of identical size.
//! \note All ids in the input ids container must have a key in the map.
//!   Otherwise an exception is thrown.
//! \note It is okay to call this function with the ids container empty but not
//!   okay to pass an empty map.
// *****************************************************************************
{
  Assert( !map.empty(), "Map must not be empty" );

  for (auto& i : ids) i = tk::cref_find( map, i );
}

std::vector< std::size_t >
remap( const std::vector< std::size_t >& ids,
       const std::unordered_map< std::size_t, std::size_t >& map )
// *****************************************************************************
//  Create remapped vector of indices using a map
//! \param[in] ids Vector of integer IDs to create new container of ids from
//! \param[in] map Hash-map of key->value creating a new order
//! \return Remapped vector of ids
//! \details This function applies a mapping (reordering) to the integer IDs
//!   passed in using the map passed in. The mapping is expressed as a hash-map
//!   of key->value pairs, where the key is the original and the value is the
//!   new ids of the mapping. The function creates and returns a new container
//!   with the remapped ids of identical size of the original ids container.
//! \note All ids in the input ids container must have a key in the map.
//!   Otherwise an exception is thrown.
//! \note It is okay to call this function with the ids container empty but not
//!   okay to pass an empty map.
// *****************************************************************************
{
  Assert( !map.empty(), "Map must not be empty" );

  // in terms of the in-place remap of a vector using a map
  auto newids = ids;
  remap( newids, map );

  return newids;
}

std::map< int, std::vector< std::size_t > >
remap( const std::map< int, std::vector< std::size_t > >& ids,
       const std::unordered_map< std::size_t, std::size_t >& map )
// *****************************************************************************
//  Create remapped map of vector of indices using a map
//! \param[in] ids Map of vector of integer IDs to create new container of ids
//!   from
//! \param[in] map Hash-map of key->value creating a new order
//! \return Remapped vector of ids
//! \details This function applies a mapping (reordering) to the map of integer
//!   IDs passed in using the map passed in by applying remap(vector,map) on
//!   each vector of ids. The keys in the returned map will be the same as in
//!   ids.
// *****************************************************************************
{
  Assert( !map.empty(), "Map must not be empty" );

  // in terms of the in-place remap of a vector using a map
  auto newids = ids;
  for (auto& m : newids) remap( m.second, map );

  return newids;
}

std::vector< std::size_t >
renumber( const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& psup )
// *****************************************************************************
//  Reorder mesh points with the advancing front technique
//! \param[in] psup Points surrounding points
//! \return Mapping created by renumbering (reordering)
// *****************************************************************************
{
  // Find out number of nodes in graph
  auto npoin = psup.second.size()-1;

  // Construct mapping using advancing front
  std::vector< int > hpoin( npoin, -1 ), lpoin( npoin, 0 );
  std::vector< std::size_t > map( npoin, 0 );
  hpoin[0] = 0;
  lpoin[0] = 1;
  std::size_t num = 1;
  while (num < npoin) {
    std::size_t cnt = 0;
    std::size_t i = 0;
    std::vector< int > kpoin( npoin, -1 );
    int p;
    while ((p = hpoin[i]) != -1) {
      ++i;
      auto P = static_cast< std::size_t >( p );
      for (auto j=psup.second[P]+1; j<=psup.second[P+1]; ++j) {
        auto q = psup.first[j];
        if (lpoin[q] != 1) {    // consider points not yet counted
          map[q] = num++;
          kpoin[cnt] = static_cast< int >( q ); // register point as counted
          lpoin[q] = 1;                         // register the point as counted
          ++cnt;
        }
      }
    }
    hpoin = kpoin;
  }

//   // Construct new->old id map
//   std::size_t i = 0;
//   std::vector< std::size_t > oldmap( npoin );
//   for (auto n : map) oldmap[n] = i++;

  // Return old->new and new->old maps
  return map;
}

std::unordered_map< std::size_t, std::size_t >
assignLid( const std::vector< std::size_t >& gid )
// *****************************************************************************
//  Assign local ids to global ids
//! \param[in] gid Global ids
//! \return Map associating global ids to local ids
// *****************************************************************************
{
  std::unordered_map< std::size_t, std::size_t > lid;
  std::size_t l = 0;
  for (auto p : gid) lid[p] = l++;
  return lid;
}

std::tuple< std::vector< std::size_t >,
            std::vector< std::size_t >,
            std::unordered_map< std::size_t, std::size_t > >
global2local( const std::vector< std::size_t >& ginpoel )
// *****************************************************************************
//  Generate element connectivity of local node IDs from connectivity of global
//  node IDs also returning the mapping between local to global IDs
//! \param[in] ginpoel Element connectivity with global node IDs
//! \return Tuple of (1) element connectivity with local node IDs, (2) the
//!   vector of unique global node IDs (i.e., the mapping between local to
//!   global node IDs), and (3) mapping between global to local node IDs.
// *****************************************************************************
{
  // Make a copy of the element connectivity with global node ids
  auto gid = ginpoel;

  // Generate a vector that holds only the unique global mesh node ids
  tk::unique( gid );

  // Assign local node ids to global node ids
  const auto lid = tk::assignLid( gid );

  Assert( gid.size() == lid.size(), "Size mismatch" );

  // Generate element connectivity using local node ids
  std::vector< std::size_t > inpoel( ginpoel.size() );
  std::size_t j = 0;
  for (auto p : ginpoel) inpoel[ j++ ] = tk::cref_find( lid, p );

  // Return element connectivty with local node IDs
  return std::make_tuple( inpoel, gid, lid );
}

bool
positiveJacobians( const std::vector< std::size_t >& inpoel,
                   const std::array< std::vector< real >, 3 >& coord )
// *****************************************************************************
// Test for positivity of the Jacobian for all cells in mesh
//! \param[in] inpoel Element connectivity (zero-based, i.e., local if parallel)
//! \param[in] coord Node coordinates
//! \return True if Jacobians of all mesh cells are positive
// *****************************************************************************
{
  Assert( !inpoel.empty(), "Mesh connectivity empty" );
  Assert( inpoel.size() % 4 == 0,
          "Mesh connectivity size must be divisible by 4 " );
  Assert( tk::uniquecopy(inpoel).size() == coord[0].size(), "Number of unique "
          "nodes in mesh connectivity must equal the number of nodes to which "
          "coordinates have been supplied" );
  Assert( tk::uniquecopy(inpoel).size() == coord[1].size(), "Number of unique "
          "nodes in mesh connectivity must equal the number of nodes to which "
          "coordinates have been supplied" );
  Assert( tk::uniquecopy(inpoel).size() == coord[2].size(), "Number of unique "
          "nodes in mesh connectivity must equal the number of nodes to which "
          "coordinates have been supplied" );
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );

  const auto& x = coord[0];
  const auto& y = coord[1];
  const auto& z = coord[2];

  for (std::size_t e=0; e<inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                           inpoel[e*4+2], inpoel[e*4+3] }};
    // compute element Jacobi determinant / (5/120) = element volume * 4
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    if (tk::triple( ba, ca, da ) < 0) return false;
 }

 return true;
}

std::map< int, std::vector< std::size_t > >
bfacenodes( const std::map< int, std::vector< std::size_t > >& bface,
            const std::vector< std::size_t >& triinpoel )
// *****************************************************************************
// Generate nodes of side set faces
//! \param[in] bface Boundary-faces mapped to side set ids
//! \param[in] triinpoel Boundary-face connectivity
//! \return Nodes of side set faces for each side set passed in
// *****************************************************************************
{
  auto bfn = bface;

  for (auto& [s,b] : bfn) {
    std::vector< std::size_t > nodes;
    for (auto f : b) {
      nodes.push_back( triinpoel[f*3+0] );
      nodes.push_back( triinpoel[f*3+1] );
      nodes.push_back( triinpoel[f*3+2] );
    }
    tk::unique( nodes );
    b = std::move( nodes );
  }

  return bfn;
}

} // tk::
