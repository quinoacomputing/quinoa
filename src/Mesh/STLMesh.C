//******************************************************************************
/*!
  \file      src/Mesh/STLMesh.C
  \author    J. Bakosi
  \date      Sun 21 Jul 2013 07:03:17 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ASCII STL (STereoLithography) mesh class definition
  \details   ASCII STL (STereoLithography) mesh class definition
*/
//******************************************************************************

#include <STLMesh.h>

using namespace Quinoa;

STLMesh::STLMesh(Memory* const memory) :
  m_memory(memory),
  m_name(""),
  m_x(),
  m_y(),
  m_z(),
  m_nodelist(),
  m_nnodes(0)
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
}

STLMesh::~STLMesh() noexcept
//******************************************************************************
//  Destructor: free memory mesh entries held and containers
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // Free memory entries held
  m_memory->freeEntry(m_nodelist);
  m_memory->freeEntry(m_z);
  m_memory->freeEntry(m_y);
  m_memory->freeEntry(m_x);
}

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

  // Allocate memory entries to store the x, y, z coordinates
  m_x = m_memory->newEntry<real>(nnodes, REAL, SCALAR, "STL x coord");
  m_y = m_memory->newEntry<real>(nnodes, REAL, SCALAR, "STL y coord");
  m_z = m_memory->newEntry<real>(nnodes, REAL, SCALAR, "STL z coord");

  // Allocate memory entry to store the node indices describing facets
  m_nodelist = m_memory->newEntry<int>(nnodes, INT, SCALAR, "STL nodelist");
  // Fill nodelist with increasing integers. This serves as the element
  // connectivity.
  for (size_t i=0; i<nnodes; ++i) m_nodelist[i] = i;
}
