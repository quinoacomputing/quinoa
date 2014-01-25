//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:26:22 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition
*/
//******************************************************************************

#include <make_unique.h>

#include <STLMesh.h>

using quinoa::STLMesh;

void
STLMesh::alloc(const size_t nnodes)
//******************************************************************************
//  Allocate memory for mesh
//! \param[in]  nnodes  Number of vertices (nodes) in mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Store number of nodes
  m_nnodes = nnodes;

  // Allocate memory to store the x, y, z coordinates
  m_x = tk::make_unique< tk::real[] >( nnodes );
  m_y = tk::make_unique< tk::real[] >( nnodes );
  m_z = tk::make_unique< tk::real[] >( nnodes );

  // Allocate memory to store the node indices describing facets
  m_nodelist = tk::make_unique< int[] >( nnodes );
  // Fill nodelist with increasing integers; this serves as connectivity
  for (size_t i=0; i<nnodes; ++i) m_nodelist[i] = i;
}
