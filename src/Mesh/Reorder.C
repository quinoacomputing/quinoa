//******************************************************************************
/*!
  \file      src/Mesh/Reorder.C
  \author    J. Bakosi
  \date      Mon 20 Apr 2015 09:58:45 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh reordering routines for unstructured meshes
  \details   Mesh reordering routines for unstructured meshes.
*/
//******************************************************************************

#include <algorithm>

#include <Reorder.h>
#include <Exception.h>

namespace tk {

void
shiftToZero( std::vector< std::size_t >& inpoel )
//******************************************************************************
//  Shift node IDs to start with zero in element connectivity
//! \param[inout] inpoel Inteconnectivity of points and elements
//! \details This function implements a simple reordering of the node ids of the
//!   element connectivity in inpoel by shifting the node ids so that the
//!   smallest is zero.
//! \note It is okay to call this function with an empty container; it will
//!    simply return without throwing an exception.
//! \author J. Bakosi
//******************************************************************************
{
  if (inpoel.empty()) return;

  // find smallest node id
  auto minId = *std::min_element( begin(inpoel), end(inpoel) );

  // shift node ids to start from zero
  for (auto& n : inpoel) n -= minId;
}

void
remap( const std::vector< std::size_t >& newid,
       std::vector< std::size_t >& inpoel )
//******************************************************************************
//  Reorder mesh points given a new order, i.e., index map
//! \param[in] newid Array of indices creating a new order
//! \param[inout] inpoel Inteconnectivity of points and elements
//! \details This function implements a simple reordering of the node ids of the
//!   element connectivity in inpoel using the vector newid. Thus the vector in
//!   newid can be thought of as a mapping between the array index to its value.
//!   The function overwrites every value, n, of inpoel with newid[ n ].
//! \note It is okay to call this function with either of the containers empty;
//!   it will simply return without throwing an exception.
//! \author J. Bakosi
//******************************************************************************
{
  if (inpoel.empty() || newid.empty()) return;

  Assert( *max_element( begin(inpoel), end(inpoel) ) < newid.size(),
          "attempt to index out of node id bounds using newid" );

  // remap node ids in element connectivity
  for (auto& n : inpoel) n = newid[ n ];
}

} // tk::
