// *****************************************************************************
/*!
  \file      src/Mesh/Reorder.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
#include <array>

#include "Types.h"

namespace tk {

//! Shift node IDs to start with zero in element connectivity
std::size_t
shiftToZero( std::vector< std::size_t >& inpoel );

//! Apply new mapping to vector of indices
void
remap( std::vector< std::size_t >& id, const std::vector< std::size_t >& map );

//! Apply new mapping to vector of real numbers
void
remap( std::vector< tk::real >& r, const std::vector< std::size_t >& map );

//! Reorder mesh points with the advancing front technique
std::vector< std::size_t >
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

//! Test positivity of the Jacobian for all cells in a mesh
bool
positiveJacobians( const std::vector< std::size_t >& inpoel,
                   const std::array< std::vector< real >, 3 >& coord );

} // ::tk

#endif // Reorder_h
