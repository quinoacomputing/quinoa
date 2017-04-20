// *****************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
// *****************************************************************************

#include <iterator>
#include <iomanip>
#include <algorithm>
#include <cstddef>
#include <ostream>
#include <string>
#include <utility>

#include "Exception.h"
#include "UnsMesh.h"
#include "StrConvUtil.h"
#include "GmshMeshWriter.h"

using tk::GmshMeshWriter;

GmshMeshWriter::GmshMeshWriter( const std::string& filename,
                                GmshFileType type,
                                tk::real version,
                                int datasize ) :
  Writer( filename ), m_type( type )
// *****************************************************************************
//  Constructor: write mandatory "$MeshFormat" section
//! \param[in] filename File to open as a Gmsh file
//! \param[in] type Gmsh file type: ASCII or binary
//! \param[in] version Gmsh file version
//! \param[in] datasize Size of double precision number on machine
//! \author J. Bakosi
// *****************************************************************************
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
GmshMeshWriter::writeMesh( const UnsMesh& mesh )
// *****************************************************************************
//  Write Gmsh mesh file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Write sections
  writeNodes( mesh );
  writeElements( mesh );
}

void
GmshMeshWriter::writeNodes( const UnsMesh& mesh )
// *****************************************************************************
//  Write "$Nodes--$EndNodes" section
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << "$Nodes" << std::endl;

  // Write out number of nodes
  m_outFile << mesh.nnode() << std::endl;

  // Write node ids and coordinates: node-number x-coord y-coord z-coord
  if (isASCII()) {
    for (std::size_t i=0; i<mesh.nnode(); ++i) {
      m_outFile << i+1 << " " << std::setprecision(16)
                << mesh.x()[i] << " "
                << mesh.y()[i] << " "
                << mesh.z()[i] << std::endl;
    }
  } else {
    for (std::size_t i=0; i<mesh.nnode(); ++i) {
      // gmsh likes one-based node ids
      int I = static_cast< int >( i+1 );
      m_outFile.write(
        reinterpret_cast<const char*>(&I), sizeof(int) );
      m_outFile.write(
        reinterpret_cast<const char*>(&mesh.x()[i]), sizeof(double) );
      m_outFile.write(
        reinterpret_cast<const char*>(&mesh.y()[i]), sizeof(double) );
      m_outFile.write(
        reinterpret_cast<const char*>(&mesh.z()[i]), sizeof(double) );
    }
    m_outFile << std::endl;
  }

  m_outFile << "$EndNodes" << std::endl;
}

void
GmshMeshWriter::writeElements( const UnsMesh& mesh )
// *****************************************************************************
//  Write "$Elements--$EndElements" section
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  m_outFile << "$Elements" << std::endl;

  // Write out number of elements
  m_outFile << mesh.lininpoel().size()/2 +
               mesh.triinpoel().size()/3 +
               mesh.tetinpoel().size()/4
            << std::endl;

  // Write out line element ids and connectivity (node list)
  writeElemBlock( 2, GmshElemType::LIN, mesh.lininpoel() );

  // Write out triangle element ids and connectivity (node list)
  writeElemBlock( 3, GmshElemType::TRI, mesh.triinpoel() );

  // Write out terahedron element ids and connectivity (node list)
  writeElemBlock( 4, GmshElemType::TET, mesh.tetinpoel() );

  if (isBinary()) m_outFile << std::endl;
  m_outFile << "$EndElements" << std::endl;
}

void
GmshMeshWriter::writeElemBlock( std::size_t nnpe,
                                GmshElemType type,
                                const std::vector< std::size_t >& inpoel )
// *****************************************************************************
//  Write element block: element ids and connectivity (node list)
//! \param[in] nnpe Number of nodes per element
//! \param[in] type Element type
//! \param[in] inpoel Element connectivity (must be zero-based)
//! \author J. Bakosi
// *****************************************************************************
{
  // Return if connectivity is empty, there is no such element block in mesh
  if (inpoel.empty()) return;

  // Make sure element connectivity starts with zero
  Assert( *std::minmax_element( begin(inpoel), end(inpoel) ).first == 0,
          "node ids should start from zero" );

  // Get number of elements in mesh
  auto n = inpoel.size()/nnpe;

  // Ignore element tags
  std::vector< std::vector< int > > tg;
  tg.resize( n );
  for (auto& t : tg) t.push_back( 0 );

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
    int nel = static_cast< int >( n );
    // elm-type num-of-elm-follow number-of-tags
    m_outFile.write( reinterpret_cast<char*>(&type), sizeof(int) );
    m_outFile.write( reinterpret_cast<char*>(&nel), sizeof(int) );
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
