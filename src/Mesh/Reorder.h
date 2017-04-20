// *****************************************************************************
/*!
  \file      src/Mesh/Reorder.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mesh reordering routines for unstructured meshes
  \details   Mesh reordering routines for unstructured meshes.
*/
// *****************************************************************************
#ifndef Reorder_h
#define Reorder_h

#include <vector>
#include <utility>
#include <tuple>
#include <unordered_map>
#include <cstddef>

namespace tk {

//! Shift node IDs to start with zero in element connectivity
std::size_t
shiftToZero( std::vector< std::size_t >& inpoel );

//! Reorder mesh point ids in a vector given a new order, i.e., index map
void
remap( std::vector< std::size_t >& id, const std::vector< std::size_t >& newid );

//! Reorder mesh points with the advancing front technique
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
renumber( const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& psup );

//! Assign local ids to global ids
std::unordered_map< std::size_t, std::size_t >
assignLid( const std::vector< std::size_t >& gid );

//! \brief Generate element connectivity of local node IDs from connectivity of
//!   global node IDs also returning the mapping between local to global IDs
std::tuple< std::vector< std::size_t >,
            std::vector< std::size_t >,
            std::unordered_map< std::size_t, std::size_t > >
global2local( const std::vector< std::size_t >& ginpoel );

} // ::tk

#endif // Reorder_h
