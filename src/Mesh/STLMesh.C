// *****************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition.
*/
// *****************************************************************************

#include <memory>

#include "STLMesh.h"

using tk::STLMesh;

void
STLMesh::alloc( std::size_t num )
// *****************************************************************************
//  Allocate memory for mesh
//! \param[in] num Number of vertices (nodes) in mesh
// *****************************************************************************
{
  // Store number of nodes
  m_nnode = num;

  // Allocate memory to store the x, y, z coordinates
  m_x = std::make_unique< tk::real[] >( num );
  m_y = std::make_unique< tk::real[] >( num );
  m_z = std::make_unique< tk::real[] >( num );
}
