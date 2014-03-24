//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.C
  \author    J. Bakosi
  \date      Mon Mar 24 15:08:27 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <limits>
#include <cmath>

#include <GmshMesh.h>
#include <GmshTxtMeshReader.h>
#include <Exception.h>

using quinoa::GmshTxtMeshReader;

GmshTxtMeshReader::GmshTxtMeshReader(const std::string filename,
                                     GmshMesh& mesh) :
  Reader(filename),
  m_mesh(mesh)
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Gmsh element types and their number of nodes,
  // all Gmsh-supported listed, Quinoa-supported at this time uncommented
  m_GmshElemNodes.insert(std::make_pair(1, 2));  // 2-node line
  m_GmshElemNodes.insert(std::make_pair(2, 3));  // 3-node triangle
  //           { 3,  4},  //! 4-node quadrangle
  //           { 4,  4},  //! 4-node tetrahedron
  //           { 5,  8},  //! 8-node hexahedron
  //           { 6,  6},  //! 6-node prism
  //           { 7,  5},  //! 5-node pyramid
  //           { 8,  3},  //! 3-node second order line
  //           { 9,  6},  //! 6-node second order triangle
  //           {10,  9},  //! 9-node second order quadrangle
  //           {11, 10},  //! 10-node second order tetrahedron
  //           {12, 27},  //! 27-node second order hexahedron
  //           {13, 18},  //! 18-node second order prism
  //           {14, 14},  //! 14-node second order pyramid
  //           {15,  1},  //! 1-node point
  //           {16,  8},  //! 8-node second order quadrangle
  //           {17, 20},  //! 20-node second order hexahedron
  //           {18, 15},  //! 15-node second order prism
  //           {19, 13},  //! 13-node second order pyramid
  //           {20,  9},  //! 9-node third order incomplete triangle
  //           {21, 10},  //! 10-node third order triangle
  //           {22, 12},  //! 12-node fourth order incomplete triangle
  //           {23, 15},  //! 15-node fourth order triangle
  //           {24, 15},  //! 15-node fifth order incomplete triangle
  //           {25, 21},  //! 21-node fifth order complete triangle
  //           {26,  4},  //! 4-node third order edge
  //           {27,  5},  //! 5-node fourth order edge
  //           {28,  6},  //! 6-node fifth order edge
  //           {29, 20},  //! 20-node third order tetrahedron
  //           {30, 35},  //! 35-node fourth order tetrahedron
  //           {31, 56},  //! 56-node fifth order tetrahedron
  //           {92, 64},  //! 64-node third order hexahedron
  //           {93,125}   //! 125-node fourth order hexahedron
}

void
GmshTxtMeshReader::read()
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
    if (s=="$Nodes") readNodes();
    else if (s=="$Elements") readElements();
    else if (s=="$PhysicalNames") readPhysicalNames();
  }

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();
}

void
GmshTxtMeshReader::readMeshFormat()
//******************************************************************************
//  Read mandatory "$MeshFormat--$EndMeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  std::string s;

  // Read in beginning of header: $MeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$MeshFormat", tk::ExceptType::FATAL,
          "Unsupported mesh format: " + m_filename );

  // Read in "version-number file-type data-size"
  tk::real version;
  int type, datasize;
  m_inFile >> version >> type >> datasize;
  ErrChk( ( fabs(version-2.2) < std::numeric_limits<tk::real>::epsilon() ||
            fabs(version-2.0) < std::numeric_limits<tk::real>::epsilon() ) &&
            type == 0 && datasize == sizeof(tk::real),
            tk::ExceptType::FATAL, "Unsupported mesh format: " + m_filename );
  getline(m_inFile, s);  // finish reading the line
  // Save version, type, datasize
  m_mesh.setVersion( version );
  m_mesh.setType( type );
  m_mesh.setDatasize( datasize );

  // Read in end of header: $EndMeshFormat
  getline( m_inFile, s );
  ErrChk( s == "$EndMeshFormat", tk::ExceptType::FATAL,
          "Unsupported mesh format: " + m_filename );
}

void
GmshTxtMeshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  std::size_t nnode;
  m_inFile >> nnode;
  ErrChk( nnode > 0, tk::ExceptType::FATAL,
          "Number of nodes must be greater than zero" );

  // Read in node ids and coordinates: node-number x-coord y-coord z-coord
  for ( std::size_t i=0; i<nnode; ++i ) {
    tk::point coord;
    int id;
    m_inFile >> id >> coord[0] >> coord[1] >> coord[2];
    m_mesh.nodeId().push_back( id );
    m_mesh.coord().push_back( coord );
  }
  std::string s;
  getline( m_inFile, s );  // finish reading the last line

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndNodes", tk::ExceptType::FATAL,
          "Unsupported mesh format: " + m_filename );
}

void
GmshTxtMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  int nel;
  m_inFile >> nel;
  ErrChk( nel > 0, tk::ExceptType::FATAL,
          "Number of elements must be greater than zero" );

  // Read in element ids, tags, and element connectivity (node list)
  for (int i=0; i<nel; ++i) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    int id, type, ntags;
    m_inFile >> id >> type >> ntags;

    // Find element type, throw exception if not supported
    auto it = m_GmshElemNodes.find( type );
    ErrChk( it != m_GmshElemNodes.end(), tk::ExceptType::FATAL,
            "Unsupported element type in mesh file: " + m_filename );

    // Read and add element tags
    std::vector<int> tags( ntags, 0 );
    for (int j=0; j<ntags; j++) {
      m_inFile >> tags[j];
    }
    addElemTags( type, tags );

    // Read and add element node list (i.e. connectivity)
    int nnode = it->second;
    std::vector< int > nodes( nnode, 0 );
    for (int j=0; j<nnode; j++) {
      m_inFile >> nodes[j];
    }
    addElem(type, nodes);

    // Put in elemId and increase counter for different types of elements
    switch (type) {
      case 1: m_mesh.lineId().push_back( id ); break;
      case 2: m_mesh.triangleId().push_back( id ); break;
    }
  }
  std::string s;
  getline( m_inFile, s );  // finish reading the last line

  // Read in end of header: $EndNodes
  getline( m_inFile, s );
  ErrChk( s == "$EndElements", tk::ExceptType::FATAL,
          "Unsupported mesh format: " + m_filename );
}

void
GmshTxtMeshReader::countPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section and count physicals
//! \author J. Bakosi
//******************************************************************************
{
  Throw(tk::ExceptType::WARNING,
      "Mesh section '$PhysicalNames -- $EndPhysicalNames' not yet implemented");
}

void
GmshTxtMeshReader::readPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
//******************************************************************************
{
  Throw(tk::ExceptType::WARNING,
      "Mesh section '$PhysicalNames -- $EndPhysicalNames' not yet implemented");
}

void
GmshTxtMeshReader::addElem( int type, const std::vector<int>& nodes )
//******************************************************************************
//  Add new element (connectivity) to given element type container
//! \param[in]  type    Elem type (container) to add to (lines, triangles, etc)
//! \param[in]  nodes   Connectivity (node Ids) of the new element
//! \author J. Bakosi
//******************************************************************************
{
  switch ( type ) {
    case 1: m_mesh.linpoel().push_back( nodes ); break;
    case 2: m_mesh.tinpoel().push_back( nodes ); break;
    default:
      Throw(tk::ExceptType::FATAL,
            "Unsupported element type in mesh file: " + m_filename);
  }
}

void
GmshTxtMeshReader::addElemTags( int type, const std::vector< int >& tags )
//******************************************************************************
//  Add new element tags to given element type
//! \param[in]  type    Elem type (container) to add to (lines, triangles, etc)
//! \param[in]  tags    Vector of tags to be added
//! \author J. Bakosi
//******************************************************************************
{
  switch ( type ) {
    case 1: m_mesh.lintag().push_back( tags ); break;
    case 2: m_mesh.tritag().push_back( tags ); break;
    default:
      Throw(tk::ExceptType::FATAL,
            "Unsupported element type in mesh file: " + m_filename);
  }
}
