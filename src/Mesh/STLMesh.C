//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Sat 30 Apr 2016 06:21:18 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition.
*/
//******************************************************************************

#include "STLMesh.h"

using tk::STLMesh;

void
STLMesh::alloc( std::size_t num )
//******************************************************************************
//  Allocate memory for mesh
//! \param[in] num Number of vertices (nodes) in mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Store number of nodes
  m_nnode = num;

  // Allocate memory to store the x, y, z coordinates
  m_x = std::make_unique< tk::real[] >( num );
  m_y = std::make_unique< tk::real[] >( num );
  m_z = std::make_unique< tk::real[] >( num );
}
