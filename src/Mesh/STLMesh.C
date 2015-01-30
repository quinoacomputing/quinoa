//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 01:00:42 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition.
*/
//******************************************************************************

#include <make_unique.h>

#include <STLMesh.h>

using quinoa::STLMesh;

void
STLMesh::alloc( const size_t nnode )
//******************************************************************************
//  Allocate memory for mesh
//! \param[in] nnode Number of vertices (nodes) in mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Store number of nodes
  m_nnodes = nnode;

  // Allocate memory to store the x, y, z coordinates
  m_x = tk::make_unique< tk::real[] >( nnode );
  m_y = tk::make_unique< tk::real[] >( nnode );
  m_z = tk::make_unique< tk::real[] >( nnode );

  // Allocate memory to store the node indices describing facets
  m_nodelist = tk::make_unique< int[] >( nnode );
  // Fill nodelist with increasing integers; this serves as connectivity
  for (size_t i=0; i<nnode; ++i) m_nodelist[i] = i;
}
