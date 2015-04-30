//******************************************************************************
/*!
  \file      src/Mesh/Reorder.h
  \author    J. Bakosi
  \date      Mon 20 Apr 2015 09:58:26 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh reordering routines for unstructured meshes
  \details   Mesh reordering routines for unstructured meshes.
*/
//******************************************************************************
#ifndef Reorder_h
#define Reorder_h

#include <vector>

namespace tk {

//! Shift node IDs to start with zero in element connectivity
void shiftToZero( std::vector< std::size_t >& inpoel );

//! Reorder mesh points given a new order, i.e., index map
void remap( const std::vector< std::size_t >& newid,
            std::vector< std::size_t >& inpoel );

} // tk::

#endif // Reorder_h
