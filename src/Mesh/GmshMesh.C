//******************************************************************************
/*!
  \file      src/Mesh/GmshMesh.C
  \author    J. Bakosi
  \date      Sat 05 Apr 2014 02:06:53 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh class definition
  \details   Gmsh mesh class definition
*/
//******************************************************************************

#include <iostream>
#include <iterator>

#include <GmshMesh.h>

using quinoa::GmshMesh;

void
GmshMesh::echoElemSets() const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  using ST = std::vector< std::vector< int > >::size_type;

  // Echo all lines
  std::cout << "* Lines elements: " << std::endl;
  std::cout << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  ST num = m_lininpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << m_lineId[i] << " " << 1 << " {";

    copy( m_lintag[i].begin(), m_lintag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_lintag[i].back() << "} {";

    copy( m_lininpoel[i].begin(), m_lininpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_lininpoel[i].back() << "}" << std::endl;
  }

  // Echo all triangles
  std::cout << "* Triangles: " << std::endl;
  std::cout << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  num = m_triinpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << m_triangleId[i] << " " << 2 << " {";

    copy( m_tritag[i].begin(), m_tritag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tritag[i].back() << "} {";

    copy( m_triinpoel[i].begin(), m_triinpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_triinpoel[i].back() << "}" << std::endl;
  }

  // Echo all tetrahedra
  std::cout << "* Tetrahedra: " << std::endl;
  std::cout << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  num = m_tetinpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << m_tetrahedronId[i] << " " << 4 << " {";

    copy( m_tettag[i].begin(), m_tettag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tettag[i].back() << "} {";

    copy( m_tetinpoel[i].begin(), m_tetinpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tetinpoel[i].back() << "}" << std::endl;
  }
}
