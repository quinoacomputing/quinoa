// *****************************************************************************
/*!
  \file      src/IO/GmshMeshReader.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
// *****************************************************************************

#include <limits>
#include <array>
#include <cmath>
#include <cstddef>
#include <istream>
#include <string>
#include <utility>
#include <vector>
#include <iostream>     // NOT NEEDED

#include "Endian.h"
#include "UnsMesh.h"
#include "GmshMeshReader.h"
#include "Reorder.h"
#include "StrConvUtil.h"

using tk::GmshMeshReader;

void
GmshMeshReader::readMesh( UnsMesh& mesh )
// *****************************************************************************
//  Public interface for read a Gmsh mesh from file
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Read in mandatory "$MeshFormat" section
  readMeshFormat();

  // Keep reading in sections until end of file. These sections can be in
  // arbitrary order, hence a while loop.
  while ( !m_inFile.eof() ) {
    std::string s;
    getline( m_inFile, s );
    if ( s == "$Nodes" )
      readNodes( mesh );
    else if ( s == "$Elements" )
      readElements( mesh );
    else if ( s == "$PhysicalNames" )
      readPhysicalNames();
  }
}

void
GmshMeshReader::readMeshFormat()
// *****************************************************************************
//  Read mandatory "$MeshFormat--$EndMeshFormat" section
//! \author J. Bakosi
// *****************************************************************************
{
  using tk::operator<<;
  std::string s;

  // Read in beginning of header: $MeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$MeshFormat",
          std::string("Unsupported mesh format '") + s + "' in file " +
          m_filename );

  // Read in "version-number file-type data-size"
  int type;
  m_inFile >> m_version >> type >> m_datasize;
  if (type == 0 )
    m_type = GmshFileType::ASCII;
  else if (type == 1 )
    m_type = GmshFileType::BINARY;

  ErrChk( ( fabs(m_version-2.2) < std::numeric_limits<tk::real>::epsilon() ||
            fabs(m_version-2.0) < std::numeric_limits<tk::real>::epsilon() ),
            std::string("Unsupported mesh version '") << m_version <<
            "' in file " << m_filename );

  ErrChk( ( m_type == GmshFileType::ASCII || m_type == GmshFileType::BINARY ),
            std::string("Unsupported mesh type '") << m_type << "' in file " <<
            m_filename );

  ErrChk( m_datasize == sizeof(tk::real),
          std::string("Unsupported mesh datasize '") << m_datasize <<
          "' in file " << m_filename );

  getline( m_inFile, s );  // finish reading the line

  // if file is binary, binary-read in binary "one"
  if ( isBinary() ) {
    int one;
    m_inFile.read( reinterpret_cast<char*>(&one), sizeof(int) );
    #ifdef __bg__
    one = tk::swap_endian< int >( one );
    #endif
    ErrChk( one == 1, "Endianness does not match in file " + m_filename );
    getline( m_inFile, s );  // finish reading the line
  }

  // Read in end of header: $EndMeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$EndMeshFormat",
          "'$EndMeshFormat' keyword is missing in file " + m_filename );
}

void
GmshMeshReader::readNodes( UnsMesh& mesh )
// *****************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  // Read in number of nodes in this node set
  std::size_t nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0,
          "Number of nodes must be greater than zero in file " + m_filename  );
  std::string s;
  if (isBinary()) getline( m_inFile, s );  // finish reading the line

  // Read in node ids and coordinates: node-number x-coord y-coord z-coord
  for ( std::size_t i=0; i<nnode; ++i ) {
    int id;
    std::array< tk::real, 3 > coord;

    if (isASCII()) {
      m_inFile >> id >> coord[0] >> coord[1] >> coord[2];
    } else {
      m_inFile.read( reinterpret_cast<char*>(&id), sizeof(int) );
      #ifdef __bg__
      id = tk::swap_endian< int >( id );
      #endif
      m_inFile.read( reinterpret_cast<char*>(coord.data()), 3*sizeof(double) );
      #ifdef __bg__
      coord[0] = tk::swap_endian< double >( coord[0] );
      coord[1] = tk::swap_endian< double >( coord[1] );
      coord[2] = tk::swap_endian< double >( coord[2] );
      #endif
    }

    mesh.x().push_back( coord[0] );
    mesh.y().push_back( coord[1] );
    mesh.z().push_back( coord[2] );
  }
  getline( m_inFile, s );  // finish reading the last line

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndNodes",
          "'$EndNodes' keyword is missing in file" + m_filename );
}

void
GmshMeshReader::readElements( UnsMesh& mesh )
// *****************************************************************************
//  Read "$Elements--$EndElements" section
//! \param[in] mesh Unstructured mesh object
//! \author J. Bakosi
// *****************************************************************************
{
  using tk::operator<<;

  std::string s;

  // Read in number of elements in this element set
  int nel;
  m_inFile >> nel;
  ErrChk( nel > 0, "Number of elements must be greater than zero in file " +
          m_filename );
  getline( m_inFile, s );  // finish reading the last line

  // Read in element ids, tags, and element connectivity (node list)
  int n=1;
  for (int i=0; i<nel; i+=n) {
    int id, elmtype, ntags;

    if (isASCII()) {
      // elm-number elm-type number-of-tags < tag > ... node-number-list
      m_inFile >> id >> elmtype >> ntags;
    } else {
      // elm-type num-of-elm-follow number-of-tags
      m_inFile.read( reinterpret_cast<char*>(&elmtype), sizeof(int) );
      m_inFile.read( reinterpret_cast<char*>(&n), sizeof(int) );
      m_inFile.read( reinterpret_cast<char*>(&ntags), sizeof(int) );
      #ifdef __bg__
      elmtype = tk::swap_endian< int >( elmtype );
      n = tk::swap_endian< int >( n );
      ntags = tk::swap_endian< int >( ntags );
      #endif
    }

    // Find element type, throw exception if not supported
    const auto it = m_elemNodes.find( elmtype );
    ErrChk( it != m_elemNodes.end(),
            std::string("Unsupported element type ") << elmtype <<
            " in mesh file: " << m_filename );

    for (int e=0; e<n; ++e) {
      // Read element id if binary
      if (isBinary()) {
        m_inFile.read( reinterpret_cast<char*>(&id), sizeof(int) );
        #ifdef __bg__
        id = tk::swap_endian< int >( id );
        #endif
      }

      // Read and ignore element tags
      std::vector< int > tags( static_cast<std::size_t>(ntags), 0 );
      if (isASCII()) {
        for (std::size_t j=0; j<static_cast<std::size_t>(ntags); j++)
          m_inFile >> tags[j];
      } else {
        m_inFile.read(
          reinterpret_cast<char*>(tags.data()),
          static_cast<std::streamsize>(
            static_cast<std::size_t>(ntags) * sizeof(int) ) );
        #ifdef __bg__
        for (auto& t : tags) t = tk::swap_endian< int >( t );
        #endif
      }

      // Read and add element node list (i.e. connectivity)
      std::size_t nnode = static_cast< std::size_t >( it->second );
      std::vector< std::size_t > nodes( nnode, 0 );
      if (isASCII()) {
        for (std::size_t j=0; j<nnode; j++)
          m_inFile >> nodes[j];
      } else {
        std::vector< int > nds( nnode, 0 );
        m_inFile.read(
          reinterpret_cast< char* >( nds.data() ),
          static_cast< std::streamsize >( nnode * sizeof(int) ) );
        #ifdef __bg__
        for (auto& j : nds) j = tk::swap_endian< int >( j );
        #endif
        for (std::size_t j=0; j<nnode; j++)
          nodes[j] = static_cast< std::size_t >( nds[j] );
      }
      // Put in element connectivity for different types of elements
      switch ( elmtype ) {
        case GmshElemType::LIN:
          for (const auto& j : nodes) mesh.lininpoel().push_back( j );
          break;
        case GmshElemType::TRI:
          for (const auto& j : nodes) mesh.triinpoel().push_back( j );
          break;
        case GmshElemType::TET:
          for (const auto& j : nodes) mesh.tetinpoel().push_back( j );
          break;
        case GmshElemType::PNT:
          break;     // ignore 1-node 'point element' type
        default: Throw( std::string("Unsupported element type ") << elmtype <<
                        " in mesh file: " << m_filename );
      }
    }
  }
  getline( m_inFile, s );  // finish reading the last line

  // Shift node IDs to start from zero (gmsh likes one-based node ids)
  shiftToZero( mesh.lininpoel() );
  shiftToZero( mesh.triinpoel() );
  shiftToZero( mesh.tetinpoel() );

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndElements",
          "'$EndElements' keyword is missing in file" + m_filename );
}

void
GmshMeshReader::readPhysicalNames()
// *****************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
// *****************************************************************************
{
  Throw( "Mesh section '$PhysicalNames -- $EndPhysicalNames' not implemented" );
}
