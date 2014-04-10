//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Thu 10 Apr 2014 09:34:36 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <iterator>
#include <iomanip>

#include <GmshMeshWriter.h>
#include <Exception.h>

using quinoa::GmshMeshWriter;

GmshMeshWriter::GmshMeshWriter( const std::string& filename,
                                GmshMesh& mesh,
                                tk::real version,
                                GmshFileType type,
                                int datasize ) :
  Writer( filename ), m_mesh( mesh ), m_type( type )
//******************************************************************************
//  Constructor: write mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  using tk::operator<<;

  // Write beginning of header: $MeshFormat
  m_outFile << "$MeshFormat\n";
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );

  // Write "version-number file-type data-size"
  m_outFile << version << " " << type << " " << datasize << "\n";
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );

  if (isBinary()) {
    int one = 1;
    m_outFile.write( reinterpret_cast<char*>(&one), sizeof(int) );
    m_outFile << std::endl;
  }

  // Write end of header: $EndMeshFormat
  m_outFile << "$EndMeshFormat" << std::endl;
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );
}

void
GmshMeshWriter::write()
//******************************************************************************
//  Write Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  // Write sections
  writeNodes();
  writeElements();
}

void
GmshMeshWriter::writeNodes()
//******************************************************************************
//  Write "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "$Nodes" << std::endl;

  // Write out number of nodes
  m_outFile << m_mesh.nnode() << std::endl;

  // Write node ids and coordinates: node-number x-coord y-coord z-coord
  if (isASCII()) {
    for (std::size_t i=0; i<m_mesh.nnode(); ++i) {
      m_outFile << m_mesh.nodeId()[i] << " " << std::setprecision(16)
                << m_mesh.coord()[i][0] << " "
                << m_mesh.coord()[i][1] << " "
                << m_mesh.coord()[i][2] << std::endl;
    }
  } else {
    for (std::size_t i=0; i<m_mesh.nnode(); ++i) {
      m_outFile.write(
        reinterpret_cast<char*>(&m_mesh.nodeId()[i]), sizeof(int) );
      m_outFile.write(
        reinterpret_cast<char*>(m_mesh.coord()[i].data()), 3*sizeof(double) );
    }
    m_outFile << std::endl;
  }

  m_outFile << "$EndNodes" << std::endl;
}

void
GmshMeshWriter::writeElements()
//******************************************************************************
//  Write "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "$Elements" << std::endl;

  // Write out number of elements
  m_outFile << m_mesh.lininpoel().size() + m_mesh.triinpoel().size() +
               m_mesh.tetinpoel().size() << std::endl;

  // Write out line element ids, tags, and connectivity (node list)
  writeElemBlock( GmshElemType::LIN, m_mesh.lineId(), m_mesh.lintag(),
                  m_mesh.lininpoel() );

  // Write out triangle element ids, tags, and connectivity (node list)
  writeElemBlock( GmshElemType::TRI, m_mesh.triangleId(), m_mesh.tritag(),
                  m_mesh.triinpoel() );

  // Write out terahedron element ids, tags, and connectivity (node list)
  writeElemBlock( GmshElemType::TET, m_mesh.tetrahedronId(), m_mesh.tettag(),
                  m_mesh.tetinpoel() );

  if (isBinary()) m_outFile << std::endl;
  m_outFile << "$EndElements" << std::endl;
}
