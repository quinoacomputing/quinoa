//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:00:23 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
                                UnsMesh& mesh,
                                GmshFileType type,
                                tk::real version,
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
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  // Write "version-number file-type data-size"
  m_outFile << version << " " << type << " " << datasize << "\n";
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );

  if (isBinary()) {
    int one = 1;
    m_outFile.write( reinterpret_cast<char*>(&one), sizeof(int) );
    m_outFile << std::endl;
  }

  // Write end of header: $EndMeshFormat
  m_outFile << "$EndMeshFormat" << std::endl;
  ErrChk( !m_outFile.bad(), "Failed to write to file: " + m_filename );
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
                << m_mesh.x()[i] << " "
                << m_mesh.y()[i] << " "
                << m_mesh.z()[i] << std::endl;
    }
  } else {
    for (std::size_t i=0; i<m_mesh.nnode(); ++i) {
      m_outFile.write(
        reinterpret_cast<char*>(&m_mesh.nodeId()[i]), sizeof(int) );
      m_outFile.write(reinterpret_cast<char*>(&m_mesh.x()[i]), sizeof(double));
      m_outFile.write(reinterpret_cast<char*>(&m_mesh.y()[i]), sizeof(double));
      m_outFile.write(reinterpret_cast<char*>(&m_mesh.z()[i]), sizeof(double));
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
  writeElemBlock( GmshElemType::LIN, m_mesh.linId(), m_mesh.lintag(),
                  m_mesh.lininpoel() );

  // Write out triangle element ids, tags, and connectivity (node list)
  writeElemBlock( GmshElemType::TRI, m_mesh.triId(), m_mesh.tritag(),
                  m_mesh.triinpoel() );

  // Write out terahedron element ids, tags, and connectivity (node list)
  writeElemBlock( GmshElemType::TET, m_mesh.tetId(), m_mesh.tettag(),
                  m_mesh.tetinpoel() );

  if (isBinary()) m_outFile << std::endl;
  m_outFile << "$EndElements" << std::endl;
}

void
GmshMeshWriter::writeElemBlock( GmshElemType type, std::vector< int >& id,
                                std::vector< std::vector< int > >& tag,
                                std::vector< std::vector< int > >& inpoel )
//******************************************************************************
//  Write element block: element ids, tags, and connectivity (node list)
//! \author J. Bakosi
//******************************************************************************
{
  if (tag.size() == 0 || id.size() == 0 || inpoel.size() == 0) return;

  auto n = inpoel.size();
  if (isASCII()) {
    for (std::size_t i=0; i<n; i++) {
      // elm-number elm-type number-of-tags < tag > ... node-number-list
      m_outFile << id[i] << " " << type << " " << tag[i].size() << " ";

      copy( tag[i].begin(), tag[i].end()-1,
            std::ostream_iterator< int >( m_outFile, " " ) );
      m_outFile << tag[i].back() << " ";

      copy( inpoel[i].begin(), inpoel[i].end()-1,
            std::ostream_iterator< int >( m_outFile, " ") );
      m_outFile << inpoel[i].back() << std::endl;
    }
  } else {
    int ntags = tag[0].size();
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&n), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
    for (std::size_t i=0; i<n; i++) {
      // element id
      m_outFile.write(
        reinterpret_cast<char*>(&id[i]), sizeof(int) );
      // element tags
      m_outFile.write( reinterpret_cast<char*>(tag[i].data()),
                       tag[i].size()*sizeof(int) );
      // element node list (i.e. connectivity)
      m_outFile.write( reinterpret_cast<char*>(inpoel[i].data()),
                       inpoel[i].size()*sizeof(int) );
    }
  }
}
