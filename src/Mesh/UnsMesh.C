//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Fri Apr 26 15:42:58 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class definition
  \details   Unstructured mesh class definition
*/
//******************************************************************************

#include <iostream>
#include <iterator>

#include <UnsMesh.h>
#include <MeshException.h>

using namespace Quinoa;

UnsMesh::UnsMesh(Memory* const memory) : m_memory(memory)
//******************************************************************************
//  Constructor: zero memory entry pointers held
//! \author J. Bakosi
//******************************************************************************
{
}

UnsMesh::~UnsMesh() noexcept
//******************************************************************************
//  Destructor: free memory mesh entries held and containers
//! \author J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    m_memory->freeEntry(m_coord);
    m_memory->freeEntry(m_nodeId);
    m_memory->freeEntry(m_lineId);
    m_memory->freeEntry(m_triangleId);
    // No exception leaves a destructor: if any of the above calls throws and
    // exception, e.g. m_coord.id points to an unallocated entry, a
    // MemoryException is thrown, caught inside here and we only emit a warning.
    // This ensures that terminate is not called and that we finish a
    // potentially already propagating exception.
#ifndef NDEBUG
  } catch (...) { cout << "WARNING: Exception in UnsMesh destructor" << endl; }
#endif // NDEBUG

  // Free containers held
  m_linpoel.clear();
  m_tinpoel.clear();
  m_lintag.clear();
  m_tritag.clear();
}

void
UnsMesh::alloc(const int nnodes, const int nlines, const int ntriangles)
//******************************************************************************
//  Allocate memory to read mesh in
//! \author J. Bakosi
//******************************************************************************
{
  // Allocate new memory entry to store the coordinates
  m_coord = m_memory->newEntry<real>(nnodes, REAL, VECTOR, COORDS_NAME);

  // Allocate new memory entry to store the node Ids
  m_nodeId = m_memory->newEntry<int>(nnodes, INT, SCALAR, NODES_NAME);

  // Allocate new memory entry to store the line element Ids
  m_lineId = m_memory->newEntry<int>(nlines, INT, SCALAR, LINES_NAME);

  // Allocate new memory entry to store the triangle element Ids
  m_triangleId = m_memory->newEntry<int>(ntriangles,INT,SCALAR,TRIANGLES_NAME);

  // Reserve capacity to store element connectivities and tags
  reserveElem(nlines, ntriangles);
}

void
UnsMesh::reserveElem(const int nlines, const int ntriangles)
//******************************************************************************
//  Reserve capacity to store element connectivities and tags
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_linpoel.reserve(nlines);
    m_tinpoel.reserve(ntriangles);
    m_lintag.reserve(nlines);
    m_tritag.reserve(ntriangles);
  } catch (bad_alloc&) {
    Throw(MemoryException,FATAL,BAD_ALLOC);
  }
}

void
UnsMesh::echoElemSets() const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  using ST = vector<vector<int>>::size_type;

  // Get pointers to the element ids
  int* linId = getLineId();
  int* triId = getTriangleId();

  // Echo all lines
  cout << "* Lines: " << endl;
  // elm-number elm-type number-of-tags < tag > ... node-number-list
  ST num = m_linpoel.size();
  for (ST i=0; i<num; i++) {
    cout << "  " << linId[i] << " " << 1 << " {";

    copy(m_lintag[i].begin(), m_lintag[i].end()-1,
         ostream_iterator<int>(cout,", "));
    cout << m_lintag[i].back() << "} {";

    copy(m_linpoel[i].begin(), m_linpoel[i].end()-1,
         ostream_iterator<int>(cout,", "));
    cout << m_linpoel[i].back() << "}" << endl;
  }

  // Echo all triangles
  cout << "* Triangles: " << endl;
  // elm-number elm-type number-of-tags < tag > ... node-number-list
  num = m_tinpoel.size();
  for (ST i=0; i<num; i++) {
    cout << "  " << triId[i] << " " << 2 << " {";

    copy(m_tritag[i].begin(), m_tritag[i].end()-1,
         ostream_iterator<int>(cout,", "));
    cout << m_tritag[i].back() << "} {";

    copy(m_tinpoel[i].begin(), m_tinpoel[i].end()-1,
         ostream_iterator<int>(cout,", "));
    cout << m_tinpoel[i].back() << "}" << endl;
  }
}
