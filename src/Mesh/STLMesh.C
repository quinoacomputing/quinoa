//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Mon Oct  7 10:20:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition
*/
//******************************************************************************

#include <STLMesh.h>

using namespace quinoa;

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
  m_x = std::unique_ptr<tk::real[]>(new tk::real [nnodes]);
  m_y = std::unique_ptr<tk::real[]>(new tk::real [nnodes]);
  m_z = std::unique_ptr<tk::real[]>(new tk::real [nnodes]);

  // Allocate memory to store the node indices describing facets
  m_nodelist = std::unique_ptr<int[]>(new int [nnodes]);
  // Fill nodelist with increasing integers; this serves as connectivity
  for (size_t i=0; i<nnodes; ++i) m_nodelist[i] = i;
}
