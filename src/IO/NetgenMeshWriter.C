//******************************************************************************
/*!
  \file      src/IO/NetgenMeshWriter.C
  \author    J. Bakosi
  \date      Sat 12 Apr 2014 07:37:37 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Netgen mesh writer
  \details   Netgen mesh writer
*/
//******************************************************************************

#include <iomanip>

#include <NetgenMeshWriter.h>

using quinoa::NetgenMeshWriter;

void
NetgenMeshWriter::write()
//******************************************************************************
//  Public interface for writing Netgen mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Write nodes
  writeNodes();
  // Write elements
  writeElements();
  // Write boundary conditions
  writeBCs();
}

void
NetgenMeshWriter::writeNodes()
//******************************************************************************
//  Write nodes
//! \author J. Bakosi
//******************************************************************************
{
  // Write iyt number of nodes
  m_outFile << m_mesh.nnode() << std::endl;

  // Write node coordinates: x-coord y-coord z-coord
  m_outFile << std::setprecision(6) << std::fixed;
  for ( std::size_t i=0; i<m_mesh.nnode(); ++i ) {
    m_outFile << '\t' << m_mesh.coord()[i][0]
              << '\t' << m_mesh.coord()[i][1]
              << '\t' << m_mesh.coord()[i][2] << std::endl;
  }
}

void
NetgenMeshWriter::writeElements()
//******************************************************************************
//  Write elements, i.e., connectivity
//! \author J. Bakosi
//******************************************************************************
{
  // Write out number of elements
  m_outFile << m_mesh.tetinpoel().size() << std::endl;

  // Write out element connectivity
  for (std::size_t i=0; i<m_mesh.tetinpoel().size(); ++i) {
    const auto& n = m_mesh.tetinpoel()[i];
    // mat n[1-4]
    m_outFile << '\t' << 1 << '\t' << n[3]+1 << '\t' << n[0]+1 << '\t'
                                   << n[1]+1 << '\t' << n[2]+1 << std::endl;
  }
}

void
NetgenMeshWriter::writeBCs()
//******************************************************************************
//  Write boundary conditions
//! \author J. Bakosi
//******************************************************************************
{
  // Write out number of boundary conditions
  m_outFile << m_mesh.bc().size() << std::endl;

  // Write out boundary conditions
  for (std::size_t i=0; i<m_mesh.bc().size(); ++i) {
    const auto& n = m_mesh.bc()[i];
    // id n[1-4]
    m_outFile << '\t' << n[0] << '\t' << n[1]+1 << '\t' << n[2]+1
              << '\t' << n[3]+1 << std::endl;
  }
}
