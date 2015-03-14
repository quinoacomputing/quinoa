//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Sat 14 Mar 2015 07:00:07 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition.
*/
//******************************************************************************

#include <make_unique.h>

#include <STLMesh.h>

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
  m_x = tk::make_unique< tk::real[] >( num );
  m_y = tk::make_unique< tk::real[] >( num );
  m_z = tk::make_unique< tk::real[] >( num );

  // Allocate memory to store the node indices describing facets
  m_nodelist = tk::make_unique< int[] >( num );
  // Fill nodelist with increasing integers; this serves as connectivity
  for (std::size_t i=0; i<num; ++i)
    m_nodelist[ i ] = static_cast< int >( i );
}
