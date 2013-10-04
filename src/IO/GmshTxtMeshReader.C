//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader.C
  \author    J. Bakosi
  \date      Thu 03 Oct 2013 08:34:18 PM MDT
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

using namespace quinoa;

GmshTxtMeshReader::GmshTxtMeshReader(const std::string filename,
                                     GmshMesh& mesh) :
  Reader(filename),
  m_mesh(mesh),
  m_GmshElemNodes(),
  m_nnodes(0),
  m_nLins(0),
  m_nTris(0),
  m_nodeCnt(0),
  m_linCnt(0),
  m_triCnt(0)
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
GmshTxtMeshReader::count()
//******************************************************************************
//  Count up elements, nodes, physicals
//! \author J. Bakosi
//******************************************************************************
{
  // Read in mandatory "$MeshFormat" section
  readMeshFormat();

  // Keep reading in sections until end of file
  while (!m_inFile.eof()) {
    std::string s;
    getline(m_inFile, s);
    if (s=="$Nodes") countNodes();
    else if (s=="$Elements") countElements();
    else if (s=="$PhysicalNames") countPhysicalNames();
  }

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inFile.clear();

  // Seek to beginning of file
  m_inFile.seekg(0, std::ios::beg);
}

void
GmshTxtMeshReader::read()
//******************************************************************************
//  Public interface for read Gmsh mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Count up elements, nodes and physical names first
  count();

  // Allocate memory to read mesh in
  m_mesh.alloc(m_nnodes, m_nLins, m_nTris);

  // Read in mandatory "$MeshFormat" section
  readMeshFormat();

  // Keep reading in sections until end of file
  while (!m_inFile.eof()) {
    std::string s;
    getline(m_inFile, s);
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
  getline(m_inFile, s);
  ErrChk(s == "$MeshFormat", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);

  // Read in "version-number file-type data-size"
  real version;
  int type, datasize;
  m_inFile >> version >> type >> datasize;
  ErrChk((fabs(version-2.2) < std::numeric_limits<real>::epsilon() ||
         fabs(version-2.0) < std::numeric_limits<real>::epsilon()) &&
         type == 0 && datasize == sizeof(real),
         ExceptType::FATAL, "Unsupported mesh format: " + m_filename);
  getline(m_inFile, s);  // finish reading the line
  // Save version, type, datasize
  m_mesh.setVersion(version);
  m_mesh.setType(type);
  m_mesh.setDatasize(datasize);

  // Read in end of header: $EndMeshFormat
  getline(m_inFile, s);
  ErrChk(s == "$EndMeshFormat", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);
}

void
GmshTxtMeshReader::countNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section and count nodes
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  int num;
  m_inFile >> num;

  // Count total number of nodes in file
  m_nnodes += num;

  // Read in node ids and coordinates throw all away
  for (int i=0; i<num; i++) {
    int n;
    real r;
    // node-number x-coord y-coord z-coord
    m_inFile >> n >> r >> r >> r;
  }
  std::string s;
  getline(m_inFile, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inFile, s);
  ErrChk(s == "$EndNodes", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);
}

void
GmshTxtMeshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  int num;
  m_inFile >> num;

  // Get pointers to node ids and coordinates
  int* nodeId = m_mesh.getNodeId();
  real* coord = m_mesh.getCoord();

  // Read in node ids and coordinates
  for (int i=0; i<num; ++i, ++m_nodeCnt) {
    // node-number x-coord y-coord z-coord
    int n3 = 3*m_nodeCnt;
    m_inFile >> nodeId[m_nodeCnt]
             >> coord[n3] >> coord[n3+1] >> coord[n3+2];
  }
  std::string s;
  getline(m_inFile, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inFile, s);
  ErrChk(s == "$EndNodes", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);
}

void
GmshTxtMeshReader::countElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section and count elements
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  int num;
  m_inFile >> num;

  // Read in element ids, tags, and element connectivity and throw all away
  for (int i=0; i<num; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    int e, type, ntags;
    m_inFile >> e >> type >> ntags;

    // Find element type, throw exception if not supported
    auto it = m_GmshElemNodes.find(type);
    ErrChk(it != m_GmshElemNodes.end(), ExceptType::FATAL,
           "Unsupported element type in mesh file: " + m_filename);

    // Read tags and throw all away
    for (int j=0; j<ntags; j++) {
      int t;
      m_inFile >> t;
    }

    // Read element node list and throw all away
    int nnodes = it->second;
    for (int j=0; j<nnodes; j++) {
      int n;
      m_inFile >> n;
    }

    // Count up different types of elements
    switch (type) {
      case 1: ++m_nLins; break;
      case 2: ++m_nTris; break;
    }
  }
  std::string s;
  getline(m_inFile, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inFile, s);
  ErrChk(s == "$EndElements", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);
}

void
GmshTxtMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  int num;
  m_inFile >> num;

  // Get pointers to the element ids
  int* linId = m_mesh.getLineId();
  int* triId = m_mesh.getTriangleId();

  // Read in element ids, tags, and element connectivity (node list)
  for (int i=0; i<num; ++i) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    int id, type, ntags;
    m_inFile >> id >> type >> ntags;

    // Find element type, throw exception if not supported
    auto it = m_GmshElemNodes.find(type);
    ErrChk(it != m_GmshElemNodes.end(), ExceptType::FATAL,
           "Unsupported element type in mesh file: " + m_filename);

    // Read and add element tags
    std::vector<int> tags(ntags,0);
    for (int j=0; j<ntags; j++) {
      m_inFile >> tags[j];
    }
    addElemTags(type, tags);

    // Read and add element node list (i.e. connectivity)
    int nnodes = it->second;
    std::vector<int> nodes(nnodes,0);
    for (int j=0; j<nnodes; j++) {
      m_inFile >> nodes[j];
    }
    addElem(type, nodes);

    // Put in elemId and increase counter for different types of elements
    switch (type) {
      case 1: linId[m_linCnt++] = id; break;
      case 2: triId[m_triCnt++] = id; break;
    }
  }
  std::string s;
  getline(m_inFile, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inFile, s);
  ErrChk(s == "$EndElements", ExceptType::FATAL,
         "Unsupported mesh format: " + m_filename);
}

void
GmshTxtMeshReader::countPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section and count physicals
//! \author J. Bakosi
//******************************************************************************
{
  Throw(ExceptType::WARNING,
      "Mesh section '$PhysicalNames -- $EndPhysicalNames' not yet implemented");
}

void
GmshTxtMeshReader::readPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
//******************************************************************************
{
  Throw(ExceptType::WARNING,
      "Mesh section '$PhysicalNames -- $EndPhysicalNames' not yet implemented");
}

void
GmshTxtMeshReader::addElem(int type, std::vector<int>& nodes)
//******************************************************************************
//  Add new element (connectivity) to given element type container
//! \param[in]  type    Elem type (container) to add to (lines, triangles, etc)
//! \param[in]  nodes   Connectivity (node Ids) of the new element
//! \author J. Bakosi
//******************************************************************************
{
  switch (type) {
    case 1: m_mesh.addLine(nodes); break;
    case 2: m_mesh.addTriangle(nodes); break;
    default:
      Throw(ExceptType::FATAL,
            "Unsupported element type in mesh file: " + m_filename);
  }
}

void
GmshTxtMeshReader::addElemTags(int type, std::vector<int>& tags)
//******************************************************************************
//  Add new element tags to given element type
//! \param[in]  type    Elem type (container) to add to (lines, triangles, etc)
//! \param[in]  tags    Vector of tags to be added
//! \author J. Bakosi
//******************************************************************************
{
  switch (type) {
    case 1: m_mesh.addLineTags(tags); break;
    case 2: m_mesh.addTriangleTags(tags); break;
    default:
      Throw(ExceptType::FATAL,
            "Unsupported element type in mesh file: " + m_filename);
  }
}
