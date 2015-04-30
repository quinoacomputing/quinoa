//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Mon 20 Apr 2015 06:15:22 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
//******************************************************************************

#include <iterator>
#include <iomanip>

#include <GmshMeshWriter.h>
#include <Exception.h>

using tk::GmshMeshWriter;

GmshMeshWriter::GmshMeshWriter( const std::string& filename,
                                const UnsMesh& mesh,
                                GmshFileType type,
                                tk::real version,
                                int datasize ) :
  Writer( filename ), m_mesh( mesh ), m_type( type )
//******************************************************************************
//  Constructor: write mandatory "$MeshFormat" section
//! \param[in] filename File to open as a Gmsh file
//! \param[in] mesh Unstructured mesh object to write data from
//! \param[in] type Gmsh file type: ASCII or binary
//! \param[in] version Gmsh file version
//! \param[in] datasize Size of double precision number on machine
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
      m_outFile << i+1 << " " << std::setprecision(16)
                << m_mesh.x()[i] << " "
                << m_mesh.y()[i] << " "
                << m_mesh.z()[i] << std::endl;
    }
  } else {
    for (std::size_t i=0; i<m_mesh.nnode(); ++i) {
      // gmsh likes one-based node ids
      int I = static_cast< int >( i+1 );
      m_outFile.write(
        reinterpret_cast<const char*>(&I), sizeof(int) );
      m_outFile.write(
        reinterpret_cast<const char*>(&m_mesh.x()[i]), sizeof(double) );
      m_outFile.write(
        reinterpret_cast<const char*>(&m_mesh.y()[i]), sizeof(double) );
      m_outFile.write(
        reinterpret_cast<const char*>(&m_mesh.z()[i]), sizeof(double) );
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
  m_outFile << m_mesh.lininpoel().size()/2 +
               m_mesh.triinpoel().size()/3 +
               m_mesh.tetinpoel().size()/4
            << std::endl;

  // Write out line element ids, tags, and connectivity (node list)
  writeElemBlock( 2, GmshElemType::LIN, m_mesh.lintag(), m_mesh.lininpoel() );

  // Write out triangle element ids, tags, and connectivity (node list)
  writeElemBlock( 3, GmshElemType::TRI, m_mesh.tritag(), m_mesh.triinpoel() );

  // Write out terahedron element ids, tags, and connectivity (node list)
  writeElemBlock( 4, GmshElemType::TET, m_mesh.tettag(), m_mesh.tetinpoel() );

  if (isBinary()) m_outFile << std::endl;
  m_outFile << "$EndElements" << std::endl;
}

void
GmshMeshWriter::writeElemBlock( std::size_t nnpe,
                                GmshElemType type,
                                const std::vector< std::vector< int > >& tag,
                                const std::vector< std::size_t >& inpoel )
//******************************************************************************
//  Write element block: element ids, tags, and connectivity (node list)
//! \param[in] nnpe Number of nodes per element
//! \param[in] type Element type
//! \param[in] tag Vectors of element tags
//! \param[in] inpoel Element connectivity (must be zero-based)
//! \author J. Bakosi
//******************************************************************************
{
  // Return if connectivity is empty, there is no such element block in mesh
  if (inpoel.empty()) return;

  // Make sure element connectivity starts with zero
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );

  // Get number of elements in mesh
  auto n = inpoel.size()/nnpe;

  // Create empty tag vector if there is no tag
  std::vector< std::vector< int > > tg;
  if (!tag.empty())
    tg = tag;
  else {
    tg.resize( n );
    for (auto& t : tg) t.push_back( 0 );
  }

  if (isASCII()) {

    for (std::size_t i=0; i<n; i++) {
      // elm-number elm-type number-of-tags < tag > ... node-number-list
      m_outFile << i+1 << " " << type << " " << tg[i].size() << " ";
      copy( tg[i].begin(), tg[i].end()-1,
            std::ostream_iterator< int >( m_outFile, " " ) );
      m_outFile << tg[i].back() << " ";

      // gmsh likes one-based node ids
      for (std::size_t k=0; k<nnpe; k++) m_outFile << inpoel[i*nnpe+k]+1 << " ";
      m_outFile << std::endl;
    }

  } else {

    int ntags = static_cast< int >( tg[0].size() );
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&n), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&ntags), sizeof(int) );
    for (std::size_t i=0; i<n; i++) {
      int I = static_cast< int >( i );
      // gmsh likes one-based node ids
      std::vector< int > Inpoel;
      for (std::size_t k=0; k<nnpe; ++k)
         Inpoel.push_back( static_cast< int >( inpoel[i*nnpe+k]+1 ) );
      // element id
      m_outFile.write( reinterpret_cast<const char*>(&I), sizeof(int) );
      // element tags
      m_outFile.write( reinterpret_cast<const char*>(tg[i].data()),
                       static_cast<std::streamsize>(tg[i].size()*sizeof(int)) );
      // element node list (i.e. connectivity)
      m_outFile.write( reinterpret_cast<const char*>(Inpoel.data()),
                       static_cast<std::streamsize>(nnpe*sizeof(int)) );
    }

  }
}
