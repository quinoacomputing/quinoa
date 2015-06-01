//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:23:28 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Netgen mesh writer class definition
  \details   Netgen mesh writer class definition. Only supports tetrahedra.
*/
//******************************************************************************

#include <iomanip>

#include "NetgenMeshWriter.h"
#include "Exception.h"

using tk::NetgenMeshWriter;

void
NetgenMeshWriter::write()
//******************************************************************************
//  Public interface for writing Netgen mesh
//! \author J. Bakosi
//******************************************************************************
{
  writeNodes();
  writeElements();
}

void
NetgenMeshWriter::writeNodes()
//******************************************************************************
//  Write nodes
//! \author J. Bakosi
//******************************************************************************
{
  // Write out number of nodes
  m_outFile << m_mesh.nnode() << std::endl;

  // Write node coordinates: x-coord y-coord z-coord
  m_outFile << std::setprecision(6) << std::fixed;
  for ( std::size_t i=0; i<m_mesh.nnode(); ++i ) {
    m_outFile << '\t' << m_mesh.x()[i]
              << '\t' << m_mesh.y()[i]
              << '\t' << m_mesh.z()[i] << std::endl;
  }
}

void
NetgenMeshWriter::writeElements()
//******************************************************************************
//  Write elements, i.e., connectivity
//! \author J. Bakosi
//******************************************************************************
{
  if (m_mesh.tetinpoel().empty()) return;

  // Make sure tetrahedron element connectivity starts with zero
  Assert( *std::minmax_element( begin(m_mesh.tetinpoel()),
                                end(m_mesh.tetinpoel()) ).first == 0,
          "tetrahedron node ids should start from zero" );

  // Get number of tetrahedra in mesh
  auto n = m_mesh.tetinpoel().size()/4;  

  // Write out number of tetrahedra
  m_outFile << n << std::endl;

  // Create empty tag vector if there is no tag
  std::vector< std::vector< int > > tg;
  if (!m_mesh.tettag().empty())
    tg = m_mesh.tettag();
  else {
    tg.resize( n );
    for (auto& t : tg) t.push_back( 0 );
  }

  // Write out tetrehadra element tags and connectivity
  for (std::size_t i=0; i<n; ++i) {
    // tag n[1-4]
    m_outFile << '\t' << tg[i][0]
              << '\t' << m_mesh.tetinpoel()[i*4+3]+1
              << '\t' << m_mesh.tetinpoel()[i*4+0]+1
              << '\t' << m_mesh.tetinpoel()[i*4+1]+1
              << '\t' << m_mesh.tetinpoel()[i*4+2]+1 << std::endl;
  }

  if (m_mesh.triinpoel().empty()) return;

  // Make sure triangle element connectivity starts with zero
  Assert( *std::minmax_element( begin(m_mesh.triinpoel()),
                                end(m_mesh.triinpoel()) ).first == 0,
          "triangle node ids should start from zero" );

  // Get number of triangles in mesh
  n = m_mesh.triinpoel().size()/3;

  // Write out number of triangles
  m_outFile << n << std::endl;

  // Create empty tag vector if there is no tag
  tg.clear();
  if (!m_mesh.tritag().empty())
    tg = m_mesh.tritag();
  else {
    tg.resize( n );
    for (auto& t : tg) t.push_back( 0 );
  }

  // Write out triangle element tags and connectivity
  for (std::size_t i=0; i<n; ++i) {
    // tag n[1-4]
    m_outFile << '\t' << tg[i][0]
              << '\t' << m_mesh.triinpoel()[i*3+0]+1
              << '\t' << m_mesh.triinpoel()[i*3+1]+1
              << '\t' << m_mesh.triinpoel()[i*3+2]+1 << std::endl;
  }
}
