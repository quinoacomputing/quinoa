//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.C
  \author    J. Bakosi
  \date      Sat 05 Apr 2014 01:57:14 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <iterator>
#include <iomanip>

#include <GmshTxtMeshWriter.h>
#include <Exception.h>

using quinoa::GmshTxtMeshWriter;

void
GmshTxtMeshWriter::write()
//******************************************************************************
//  Write Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  // Write out mandatory "$MeshFormat" section
  writeMeshFormat();

  // Write sections
  writeNodes();
  writeElements();
}

void
GmshTxtMeshWriter::writeMeshFormat()
//******************************************************************************
//  Write mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  // Write beginning of header: $MeshFormat
  m_outFile << "$MeshFormat\n";
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );

  // Write "version-number file-type data-size"
  m_outFile << m_mesh.getVersion() << " "
            << m_mesh.getType() << " "
            << m_mesh.getDatasize() << "\n";
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );

  // Write end of header: $EndMeshFormat
  m_outFile << "$EndMeshFormat" << std::endl;
  ErrChk( !m_outFile.bad(), tk::ExceptType::FATAL,
          "Failed to write to file: " + m_filename );
}

void
GmshTxtMeshWriter::writeNodes()
//******************************************************************************
//  Write "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  m_outFile << "$Nodes" << std::endl;

  // Write out number of nodes
  m_outFile << m_mesh.nnode() << std::endl;

  // Write node ids and coordinates
  for (std::size_t i=0; i<m_mesh.nnode(); ++i) {
    // node-number x-coord y-coord z-coord
    m_outFile << m_mesh.nodeId()[i] << " " << std::setprecision(16)
              << m_mesh.coord()[i][0] << " "
              << m_mesh.coord()[i][1] << " "
              << m_mesh.coord()[i][2] << std::endl;
  }

  m_outFile << "$EndNodes" << std::endl;
}

void
GmshTxtMeshWriter::writeElements()
//******************************************************************************
//  Write "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  using ST = std::vector< std::vector< int > >::size_type;

  m_outFile << "$Elements" << std::endl;

  // Get pointers to the element ids, connectivities, and tags
  auto& lininpoel = m_mesh.lininpoel();
  auto& lintag = m_mesh.lintag();
  auto& triinpoel = m_mesh.triinpoel();
  auto& tritag = m_mesh.tritag();
  auto& tetinpoel = m_mesh.tetinpoel();
  auto& tettag = m_mesh.tettag();

  // Get number of elements
  ST nlines = lininpoel.size();
  ST ntriangles = triinpoel.size();
  ST ntetrahedra = tetinpoel.size();

  // Write out number of elements
  m_outFile << nlines + ntriangles + ntetrahedra << std::endl;

  // Write out line element ids, tags, and connectivity (node list)
  for (ST i=0; i<nlines; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outFile << m_mesh.lineId()[i] << " " << 1 << " " << lintag[i].size()
              << " ";

    copy( lintag[i].begin(), lintag[i].end()-1,
          std::ostream_iterator< int >( m_outFile, " " ) );
    m_outFile << lintag[i].back() << " ";

    copy( lininpoel[i].begin(), lininpoel[i].end()-1,
          std::ostream_iterator< int >( m_outFile, " ") );
    m_outFile << lininpoel[i].back() << std::endl;
  }

  // Write out triangle element ids, tags, and connectivity (node list)
  for (ST i=0; i<ntriangles; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outFile << m_mesh.triangleId()[i] << " " << 2 << " " << tritag[i].size()
              << " ";

    copy(tritag[i].begin(), tritag[i].end()-1,
         std::ostream_iterator<int>(m_outFile," "));
    m_outFile << tritag[i].back() << " ";

    copy(triinpoel[i].begin(), triinpoel[i].end()-1,
         std::ostream_iterator<int>(m_outFile," "));
    m_outFile << triinpoel[i].back() << std::endl;
  }

  // Write out terahedron element ids, tags, and connectivity (node list)
  for (ST i=0; i<ntetrahedra; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outFile << m_mesh.tetrahedronId()[i] << " " << 4 << " "
              << tettag[i].size() << " ";

    copy(tettag[i].begin(), tettag[i].end()-1,
         std::ostream_iterator<int>(m_outFile," "));
    m_outFile << tettag[i].back() << " ";

    copy(tetinpoel[i].begin(), tetinpoel[i].end()-1,
         std::ostream_iterator<int>(m_outFile," "));
    m_outFile << tetinpoel[i].back() << std::endl;
  }

  m_outFile << "$EndElements" << std::endl;
}
