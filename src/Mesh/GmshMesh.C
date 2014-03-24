//******************************************************************************
/*!
  \file      src/Mesh/GmshMesh.C
  \author    J. Bakosi
  \date      Mon Mar 24 13:51:41 2014
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
  ST num = m_linpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << m_lineId[i] << " " << 1 << " {";

    copy( m_lintag[i].begin(), m_lintag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_lintag[i].back() << "} {";

    copy( m_linpoel[i].begin(), m_linpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_linpoel[i].back() << "}" << std::endl;
  }

  // Echo all triangles
  std::cout << "* Triangles: " << std::endl;
  std::cout << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  num = m_tinpoel.size();
  for (ST i=0; i<num; i++) {
    std::cout << "  " << m_triangleId[i] << " " << 2 << " {";

    copy( m_tritag[i].begin(), m_tritag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tritag[i].back() << "} {";

    copy( m_tinpoel[i].begin(), m_tinpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tinpoel[i].back() << "}" << std::endl;
  }
}
