//******************************************************************************
/*!
  \file      src/Mesh/GmshMesh.C
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:29:29 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh class definition
  \details   Gmsh mesh class definition
*/
//******************************************************************************

#include <iostream>
#include <iterator>

#include <make_unique.h>

#include <GmshMesh.h>
#include <Exception.h>

using quinoa::GmshMesh;

void
GmshMesh::alloc(const int nnodes, const int nlines, const int ntriangles)
//******************************************************************************
//  Allocate memory to read mesh in
//! \author J. Bakosi
//******************************************************************************
{
  // Store number of nodes
  m_nnodes = nnodes;

  // Allocate memory to store the coordinates
  m_coord = tk::make_unique< tk::real[] >( 3*nnodes );

  // Allocate memory to store the node Ids
  m_nodeId = tk::make_unique< int[] >( nnodes );

  // Allocate memory to store the line element Ids
  m_lineId = tk::make_unique< int[] >( nlines );

  // Allocate memory to store the triangle element Ids
  m_triangleId = tk::make_unique< int[] >( ntriangles );

  // Reserve capacity to store element connectivities and tags
  reserveElem(nlines, ntriangles);
}

void
GmshMesh::reserveElem(const int nlines, const int ntriangles) noexcept
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
      Throw(tk::ExceptType::FATAL, "Cannot allocate memory for mesh");
    }
    catch (...) {
      Throw(tk::ExceptType::UNCAUGHT,
            "non-standard exception in GmshMesh::reserveElem()");
    }
}

void
GmshMesh::echoElemSets() const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  using ST = std::vector< std::vector< int > >::size_type;

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
