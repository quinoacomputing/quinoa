//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:18:55 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Netgen mesh writer class definition
  \details   Netgen mesh writer class definition. Only supports tetrahedra.
*/
//******************************************************************************

#include <iomanip>

#include <NetgenMeshWriter.h>
#include <Exception.h>

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
  // Write out number of tetrahedra
  m_outFile << m_mesh.tetinpoel().size() << std::endl;

  // Write out tetrehadra element tags and connectivity
  for (std::size_t i=0; i<m_mesh.tetinpoel().size(); ++i) {
    // tag n[1-4]
    m_outFile << '\t' << m_mesh.tettag()[i][0]
              << '\t' << m_mesh.tetinpoel()[i][3]
              << '\t' << m_mesh.tetinpoel()[i][0]
              << '\t' << m_mesh.tetinpoel()[i][1]
              << '\t' << m_mesh.tetinpoel()[i][2] << std::endl;
  }

  // Write out number of triangles
  m_outFile << m_mesh.triinpoel().size() << std::endl;

  // Write out triangle element tags and connectivity
  for (std::size_t i=0; i<m_mesh.triinpoel().size(); ++i) {
    // tag n[1-4]
    m_outFile << '\t' << m_mesh.tritag()[i][0]
              << '\t' << m_mesh.triinpoel()[i][0]
              << '\t' << m_mesh.triinpoel()[i][1]
              << '\t' << m_mesh.triinpoel()[i][2] << std::endl;
  }
}
