//******************************************************************************
/*!
  \file      src/Mesh/Reorder.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 10:06:58 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh reordering routines for unstructured meshes
  \details   Mesh reordering routines for unstructured meshes.
*/
//******************************************************************************
#ifndef Reorder_h
#define Reorder_h

#include <vector>
#include <utility>
#include <cstddef>

namespace tk {

//! Shift node IDs to start with zero in element connectivity
void
shiftToZero( std::vector< std::size_t >& inpoel );

//! Reorder mesh point ids in a vector given a new order, i.e., index map
void
remap( std::vector< std::size_t >& id, const std::vector< std::size_t >& newid );

//! Reorder mesh points with the advancing front technique
std::pair< std::vector< std::size_t >, std::vector< std::size_t > >
renumber( std::pair< std::vector< std::size_t >,
                     std::vector< std::size_t > >& psup );

} // tk::

#endif // Reorder_h
