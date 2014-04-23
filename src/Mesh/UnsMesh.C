//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Wed Apr 23 13:39:33 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     3D unstructured mesh class definition
  \details   3D unstructured mesh class definition
*/
//******************************************************************************

#include <iostream>
#include <iterator>

#include <UnsMesh.h>

using quinoa::UnsMesh;

void
UnsMesh::echoElemSets() const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  // Echo all lines
  std::cout << "* Lines elements: " << std::endl;
  std::cout << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  std::size_t num = m_lininpoel.size();
  for (std::size_t i=0; i<num; i++) {
    std::cout << "  " << m_linId[i] << " " << 1 << " {";

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
  for (std::size_t i=0; i<num; i++) {
    std::cout << "  " << m_triId[i] << " " << 2 << " {";

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
  for (std::size_t i=0; i<num; i++) {
    std::cout << "  " << m_tetId[i] << " " << 4 << " {";

    copy( m_tettag[i].begin(), m_tettag[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tettag[i].back() << "} {";

    copy( m_tetinpoel[i].begin(), m_tetinpoel[i].end()-1,
          std::ostream_iterator< int >( std::cout, ", " ) );
    std::cout << m_tetinpoel[i].back() << "}" << std::endl;
  }
}
