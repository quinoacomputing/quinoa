//******************************************************************************
/*!
  \file      src/IO/GmshMeshReader.C
  \author    J. Bakosi
  \date      Mon 26 May 2014 04:43:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <limits>
#include <cmath>
#include <array>

#include <UnsMesh.h>
#include <GmshMeshReader.h>

using quinoa::GmshMeshReader;

void
GmshMeshReader::read()
//******************************************************************************
//  Public interface for read Gmsh mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Read in mandatory "$MeshFormat" section
  readMeshFormat();

  // Keep reading in sections until end of file
  while (!m_inFile.eof()) {
    std::string s;
    getline( m_inFile, s );
    if (s=="$Nodes")
      readNodes();
    else if (s=="$Elements")
      readElements();
    else if (s=="$PhysicalNames")
      readPhysicalNames();
  }

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();
}

void
GmshMeshReader::readMeshFormat()
//******************************************************************************
//  Read mandatory "$MeshFormat--$EndMeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  std::string s;

  // Read in beginning of header: $MeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$MeshFormat", tk::ExceptType::FATAL,
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
            tk::ExceptType::FATAL,
            std::string("Unsupported mesh version '") + m_version +
            "' in file " + m_filename );

  ErrChk( ( m_type == GmshFileType::ASCII || m_type == GmshFileType::BINARY ),
          tk::ExceptType::FATAL,
          std::string("Unsupported mesh type '") + m_type + "' in file " +
            m_filename );

  ErrChk( m_datasize == sizeof(tk::real),
          tk::ExceptType::FATAL,
          std::string("Unsupported mesh datasize '") + m_datasize +
          "' in file " + m_filename );

  getline( m_inFile, s );  // finish reading the line

  // if file is binary, binary-read in binary "one"
  if (isBinary()) {
    int one;
    m_inFile.read( reinterpret_cast<char*>(&one), sizeof(int) );
    ErrChk( one == 1,
            tk::ExceptType::FATAL,
            "Endianness does not match in file " + m_filename );
    getline( m_inFile, s );  // finish reading the line
  }

  // Read in end of header: $EndMeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$EndMeshFormat", tk::ExceptType::FATAL,
          "'$EndMeshFormat' keyword is missing in file " + m_filename );
}

void
GmshMeshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  std::size_t nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0, tk::ExceptType::FATAL,
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
      m_inFile.read( reinterpret_cast<char*>(coord.data()), 3*sizeof(double) );
    }

    m_mesh.nodeId().push_back( id );
    m_mesh.x().push_back( coord[0] );
    m_mesh.y().push_back( coord[1] );
    m_mesh.z().push_back( coord[2] );
  }
  getline( m_inFile, s );  // finish reading the last line

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndNodes", tk::ExceptType::FATAL,
          "'$EndNodes' keyword is missing in file" + m_filename );
}

void
GmshMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  using tk::operator+;

  std::string s;

  // Read in number of elements in this element set
  int nel;
  m_inFile >> nel;
  ErrChk( nel > 0, tk::ExceptType::FATAL,
          "Number of elements must be greater than zero in file " +
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
    }

    // Find element type, throw exception if not supported
    const auto it = m_elemNodes.find( elmtype );
    ErrChk( it != m_elemNodes.end(), tk::ExceptType::FATAL,
            std::string("Unsupported element type ") + elmtype +
            " in mesh file: " + m_filename );

    for (int e=0; e<n; ++e) {
      // Read element id if binary
      if (isBinary()) {
        m_inFile.read( reinterpret_cast<char*>(&id), sizeof(int) );
      }

      // Read and add element tags
      std::vector<int> tags( ntags, 0 );
      if (isASCII()) {
        for (int j=0; j<ntags; j++)
          m_inFile >> tags[j];
      } else {
        m_inFile.read(
          reinterpret_cast<char*>(tags.data()), ntags*sizeof(int) );
      }
      // Put in element tags for different types of elements
      push_back( elmtype, tags, m_mesh.lintag(), m_mesh.tritag(),
                 m_mesh.tettag() );

      // Read and add element node list (i.e. connectivity)
      int nnode = it->second;
      std::vector< int > nodes( nnode, 0 );
      if (isASCII()) {
        for (int j=0; j<nnode; j++)
          m_inFile >> nodes[j];
      } else {
        m_inFile.read(
          reinterpret_cast<char*>(nodes.data()), nnode*sizeof(int) );
      }
      // Put in element connectivity for different types of elements
      push_back( elmtype, nodes, m_mesh.lininpoel(), m_mesh.triinpoel(),
                 m_mesh.tetinpoel() );

      // Put in elemId for different types of elements
      push_back( elmtype, id, m_mesh.linId(), m_mesh.triId(), m_mesh.tetId() );
    }
  }
  getline( m_inFile, s );  // finish reading the last line

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndElements", tk::ExceptType::FATAL,
          "'$EndElements' keyword is missing in file" + m_filename );
}

void
GmshMeshReader::readPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
//******************************************************************************
{
  Throw(tk::ExceptType::WARNING,
      "Mesh section '$PhysicalNames -- $EndPhysicalNames' not yet implemented");
}
