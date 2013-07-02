//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Tue Jul  2 15:14:42 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Unstructured mesh class definition
  \details   Unstructured mesh class definition
*/
//******************************************************************************

#include <iostream>
#include <iterator>

#include <UnsMesh.h>

using namespace Quinoa;

UnsMesh::UnsMesh(Memory* const memory) :
  m_memory(memory),
  m_type(0),
  m_datasize(0),
  m_coord(),
  m_nodeId(),
  m_lineId(),
  m_triangleId(),
  m_linpoel(),
  m_tinpoel(),
  m_lintag(),
  m_tritag()
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
  m_memory->freeEntry(m_coord);
  m_memory->freeEntry(m_nodeId);
  m_memory->freeEntry(m_lineId);
  m_memory->freeEntry(m_triangleId);

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

  } catch (std::exception&) {
      Throw(ExceptType::FATAL, "Cannot allocate memory for mesh");
    }
    catch (...) {
      Throw(ExceptType::UNCAUGHT,
            "non-standard exception in UnsMesh::reserveElem()");
    }
}

void
UnsMesh::echoElemSets() const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  using ST = std::vector<std::vector<int>>::size_type;

  // Get pointers to the element ids
  int* linId = getLineId();
  int* triId = getTriangleId();

  // Echo all lines
  std::cout << "* Lines: " << std::endl;
  // elm-number elm-type number-of-tags < tag > ... node-number-list
  ST num = m_linpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << linId[i] << " " << 1 << " {";

    copy(m_lintag[i].begin(), m_lintag[i].end()-1,
         std::ostream_iterator<int>(std::cout,", "));
    std::cout << m_lintag[i].back() << "} {";

    copy(m_linpoel[i].begin(), m_linpoel[i].end()-1,
         std::ostream_iterator<int>(std::cout,", "));
    std::cout << m_linpoel[i].back() << "}" << std::endl;
  }

  // Echo all triangles
  std::cout << "* Triangles: " << std::endl;
  // elm-number elm-type number-of-tags < tag > ... node-number-list
  num = m_tinpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << triId[i] << " " << 2 << " {";

    copy(m_tritag[i].begin(), m_tritag[i].end()-1,
         std::ostream_iterator<int>(std::cout,", "));
    std::cout << m_tritag[i].back() << "} {";

    copy(m_tinpoel[i].begin(), m_tinpoel[i].end()-1,
         std::ostream_iterator<int>(std::cout,", "));
    std::cout << m_tinpoel[i].back() << "}" << std::endl;
  }
}
