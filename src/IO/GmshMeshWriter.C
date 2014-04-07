//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Mon 07 Apr 2014 07:14:14 AM MDT
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
                                int type, 
                                int datasize ) :
  Writer(filename), m_mesh(mesh), m_type(type)
//******************************************************************************
//  Constructor: write mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
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

  // Get pointers to the element ids, connectivities, and tags
  auto& lininpoel = m_mesh.lininpoel();
  auto& lintag = m_mesh.lintag();
  auto& triinpoel = m_mesh.triinpoel();
  auto& tritag = m_mesh.tritag();
  auto& tetinpoel = m_mesh.tetinpoel();
  auto& tettag = m_mesh.tettag();

  // Get number of elements
  std::size_t nlines = lininpoel.size();
  std::size_t ntriangles = triinpoel.size();
  std::size_t ntetrahedra = tetinpoel.size();

  // Write out number of elements
  m_outFile << nlines + ntriangles + ntetrahedra << std::endl;

  // Write out line element ids, tags, and connectivity (node list)
  if (isASCII()) {
    for (std::size_t i=0; i<nlines; i++) {
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
  } else {
    int type = 1;
    int ntags = lintag[0].size();
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&nlines), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
    for (std::size_t i=0; i<nlines; i++) {
      // element id
      m_outFile.write(
        reinterpret_cast<char*>(&m_mesh.lineId()[i]), sizeof(int) );
      // element tags
      m_outFile.write( reinterpret_cast<char*>(lintag[i].data()),
                       lintag[i].size()*sizeof(int) );
      // element node list (i.e. connectivity)
      m_outFile.write( reinterpret_cast<char*>(lininpoel[i].data()),
                       lininpoel[i].size()*sizeof(int) );
    }
  }

  // Write out triangle element ids, tags, and connectivity (node list)
  if (isASCII()) {
    for (std::size_t i=0; i<ntriangles; i++) {
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
  } else {
    int type = 2;
    int ntags = tritag[0].size();
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntriangles), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
    for (std::size_t i=0; i<ntriangles; i++) {
      // element id
      m_outFile.write(
        reinterpret_cast<char*>(&m_mesh.triangleId()[i]), sizeof(int) );
      // element tags
      m_outFile.write( reinterpret_cast<char*>(tritag[i].data()),
                       tritag[i].size()*sizeof(int) );
      // element node list (i.e. connectivity)
      m_outFile.write( reinterpret_cast<char*>(triinpoel[i].data()),
                       triinpoel[i].size()*sizeof(int) );
    }
  }

  // Write out terahedron element ids, tags, and connectivity (node list)
  if (isASCII()) {
    for (std::size_t i=0; i<ntetrahedra; i++) {
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
  } else {
    int type = 4;
    int ntags = tettag[0].size();
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntetrahedra), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
    for (std::size_t i=0; i<ntetrahedra; i++) {
      // element id
      m_outFile.write(
        reinterpret_cast<char*>(&m_mesh.tetrahedronId()[i]), sizeof(int) );
      // element tags
      m_outFile.write( reinterpret_cast<char*>(tettag[i].data()),
                       tettag[i].size()*sizeof(int) );
      // element node list (i.e. connectivity)
      m_outFile.write( reinterpret_cast<char*>(tetinpoel[i].data()),
                       tetinpoel[i].size()*sizeof(int) );
    }
  }

  if (isBinary()) m_outFile << std::endl;
  m_outFile << "$EndElements" << std::endl;
}
