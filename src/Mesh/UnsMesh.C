//******************************************************************************
/*!
  \file      src/Mesh/UnsMesh.C
  \author    J. Bakosi
  \date      Thu 29 May 2014 06:09:24 AM MDT
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
UnsMesh::echoElemSets( const tk::Print& print ) const
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  // Echo all lines
  print << "* Lines elements: " << std::endl;
  print << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  std::size_t num = m_lininpoel.size();
  for (std::size_t i=0; i<num; i++) {
    print << "  " << m_linId[i] << " " << 1 << " {";

    copy( m_lintag[i].begin(), m_lintag[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_lintag[i].back() << "} {";

    copy( m_lininpoel[i].begin(), m_lininpoel[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_lininpoel[i].back() << "}" << std::endl;
  }

  // Echo all triangles
  print << "* Triangles: " << std::endl;
  print << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  num = m_triinpoel.size();
  for (std::size_t i=0; i<num; i++) {
    print << "  " << m_triId[i] << " " << 2 << " {";

    copy( m_tritag[i].begin(), m_tritag[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_tritag[i].back() << "} {";

    copy( m_triinpoel[i].begin(), m_triinpoel[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_triinpoel[i].back() << "}" << std::endl;
  }

  // Echo all tetrahedra
  print << "* Tetrahedra: " << std::endl;
  print << "  elm-num elm-type {tags} {nodelist} " << std::endl;
  num = m_tetinpoel.size();
  for (std::size_t i=0; i<num; i++) {
    print << "  " << m_tetId[i] << " " << 4 << " {";

    copy( m_tettag[i].begin(), m_tettag[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_tettag[i].back() << "} {";

    copy( m_tetinpoel[i].begin(), m_tetinpoel[i].end()-1,
          std::ostream_iterator< int >( print, ", " ) );
    print << m_tetinpoel[i].back() << "}" << std::endl;
  }
}
