//******************************************************************************
/*!
  \file      src/Mesh/GmshTxtMeshReader.C
  \author    J. Bakosi
  \date      Sat 06 Oct 2012 08:24:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

#include <QuinoaTypes.h>
#include <GmshTxtMeshReader.h>
#include <MeshException.h>
#include <MemoryException.h>
#include <Memory.h>

using namespace Quinoa;

GmshTxtMeshReader::~GmshTxtMeshReader()
//******************************************************************************
//  Destructor: free mesh entries
//! \author J. Bakosi
//******************************************************************************
{
  // Free node sets
  if (m_meshEntry.size()) {
    // Free all mesh sets
    MeshSet::const_iterator it;
    for (it=m_meshEntry.begin(); it!=m_meshEntry.end(); it++) {
      m_memory->freeEntry(*it);
    }
    // Clear container
    m_meshEntry.clear();
  }
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
  while (!m_inMesh.eof()) {
    string s;
    getline(m_inMesh, s);
    if (s=="$Nodes") countNodes();
    else if (s=="$Elements") countElements();
    else if (s=="$PhysicalNames") countPhysicalNames();
  }

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inMesh.clear();

  // Seek to beginning of file
  m_inMesh.seekg(0, ios::beg);
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
  alloc();

  // Read in mandatory "$MeshFormat" section
  readMeshFormat();

  // Keep reading in sections until end of file
  while (!m_inMesh.eof()) {
    string s;
    getline(m_inMesh, s);
    if (s=="$Nodes") readNodes();
    else if (s=="$Elements") readElements();
    else if (s=="$PhysicalNames") readPhysicalNames();
  }

  // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
  m_inMesh.clear();

  echoElemSets();
}

void
GmshTxtMeshReader::echoElemSets()
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  // Echo all lines
  cout << "* Lines: " << endl;

  // elm-number elm-type number-of-tags < tag > ... node-number-list
  typedef vector<vector<Int>>::size_type ST;
  ST num = m_lin.size();
  for (ST i=0; i<num; i++) {
    cout << "  " << m_linId[i] << " " << 1 << " {";

    copy(m_linTag[i].begin(), m_linTag[i].end()-1,
         ostream_iterator<Int>(cout,", "));
    cout << m_linTag[i].back() << "} {";

    copy(m_lin[i].begin(), m_lin[i].end()-1,
         ostream_iterator<Int>(cout,", "));
    cout << m_lin[i].back() << "}" << endl;
  }

  // Echo all triangles
  cout << "* Triangles: " << endl;

  // elm-number elm-type number-of-tags < tag > ... node-number-list
  num = m_tri.size();
  for (ST i=0; i<num; i++) {
    cout << "  " << m_triId[i] << " " << 2 << " {";

    copy(m_triTag[i].begin(), m_triTag[i].end()-1,
         ostream_iterator<Int>(cout,", "));
    cout << m_triTag[i].back() << "} {";

    copy(m_tri[i].begin(), m_tri[i].end()-1,
         ostream_iterator<Int>(cout,", "));
    cout << m_tri[i].back() << "}" << endl;
  }
}

void
GmshTxtMeshReader::readMeshFormat()
//******************************************************************************
//  Read mandatory "$MeshFormat--$EndMeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  string s;

  // Read in beginning of header: $MeshFormat
  getline(m_inMesh, s);
  if (s!="$MeshFormat")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT, 
                        m_filename);

  // Read in "version-number file-type data-size"
  Real version;
  Int type, datasize;
  m_inMesh >> version >> type >> datasize;
  if ((version!=2.2 && version!=2.0) || type!=0 || datasize!=sizeof(Real))
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
  getline(m_inMesh, s);  // finish reading the line
  // Save version, type, datasize
  m_mesh->setVersion(version);
  m_mesh->setType(type);
  m_mesh->setDatasize(datasize);

  // Read in end of header: $EndMeshFormat
  getline(m_inMesh, s);
  if (s!="$EndMeshFormat")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
}

void
GmshTxtMeshReader::countNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section and count nodes
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  Int num;
  m_inMesh >> num;

  // Count total number of nodes in file
  m_nnodes += num;

  // Read in node ids and coordinates throw all away
  for (Int i=0; i<num; i++) {
    Int n;
    Real r;
    // node-number x-coord y-coord z-coord
    m_inMesh >> n >> r >> r >> r;
  }
  string s;
  getline(m_inMesh, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inMesh, s);
  if (s!="$EndNodes")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
}

void
GmshTxtMeshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  Int num;
  m_inMesh >> num;

  // Read in node ids and coordinates
  for (Int i=0; i<num; ++i, ++m_nodeCnt) {
    // node-number x-coord y-coord z-coord
    Int n3 = 3*m_nodeCnt;
    m_inMesh >> m_node[m_nodeCnt]
             >> m_coord[n3] >> m_coord[n3+1] >> m_coord[n3+2];
  }
  string s;
  getline(m_inMesh, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inMesh, s);
  if (s!="$EndNodes")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
}

void
GmshTxtMeshReader::countElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section and count elements
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  Int num;
  m_inMesh >> num;

  // Count total number of elements in file
  m_nelems += num;

  // Read in element ids, tags, and element connectivity and throw all away
  for (Int i=0; i<num; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    Int e, type, ntags;
    m_inMesh >> e >> type >> ntags;

    // Find element type, throw exception if not supported
    auto it = GmshElemNodes.find(type);
    if (it==GmshElemNodes.end())
      throw MeshException(ExceptType::FATAL,
            MeshExceptType::BAD_ELEMENT,
            m_filename);

    // Count up different types of elements
    switch (type) {
      case 1: ++m_nLins; break;
      case 2: ++m_nTris; break;
    }

    // Read tags and throw all away
    for (Int j=0; j<ntags; j++) {
      Int t;
      m_inMesh >> t;
    }

    // Read element node list and throw all away
    Int nnodes = it->second;
    for (Int j=0; j<nnodes; j++) {
      Int n;
      m_inMesh >> n;
    }
  }
  string s;
  getline(m_inMesh, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inMesh, s);
  if (s!="$EndElements")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
}

void
GmshTxtMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  Int num;
  m_inMesh >> num;

  // Read in element ids, tags, and element connectivity (node list)
  for (Int i=0; i<num; ++i) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    Int id, type, ntags;
    m_inMesh >> id >> type >> ntags;

    // Find element type, throw exception if not supported
    auto it = GmshElemNodes.find(type);
    if (it==GmshElemNodes.end())
      throw MeshException(ExceptType::FATAL,
            MeshExceptType::BAD_ELEMENT,
            m_filename);

    // Read tags
    vector<Int> tags(ntags,0);
    for (Int j=0; j<ntags; j++) {
      m_inMesh >> tags[j];
    }

    // Add tags for different types of elements
    switch (type) {
      case 1: addElemTags(m_linTag, tags); break;
      case 2: addElemTags(m_triTag, tags); break;
    }

    // Read element node list
    Int nnodes = it->second;
    vector<Int> nodes(nnodes,0);
    for (Int j=0; j<nnodes; j++) {
      m_inMesh >> nodes[j];
    }

    // Add node list for different types of elements
    switch (type) {
      case 1: addElem(m_lin, nodes); break;
      case 2: addElem(m_tri, nodes); break;
    }

    // Put in elemId and increase counter for different types of elements
    switch (type) {
      case 1: m_linId[m_linCnt++] = id; break;
      case 2: m_triId[m_triCnt++] = id; break;
    }
  }
  string s;
  getline(m_inMesh, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inMesh, s);
  if (s!="$EndElements")
    throw MeshException(ExceptType::FATAL,
                        MeshExceptType::BAD_FORMAT,
                        m_filename);
}

void
GmshTxtMeshReader::countPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section and count physicals
//! \author J. Bakosi
//******************************************************************************
{
  throw MeshException(ExceptType::WARNING,
                      MeshExceptType::UNIMPLEMENTED,
                      "$PhysicalNames--$EndPhysicalNames");
}

void
GmshTxtMeshReader::readPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
//******************************************************************************
{
  throw MeshException(ExceptType::WARNING,
                      MeshExceptType::UNIMPLEMENTED,
                      "$PhysicalNames--$EndPhysicalNames");
}

void
GmshTxtMeshReader::alloc()
//******************************************************************************
//  Allocate memory to read mesh in
//! \author J. Bakosi
//******************************************************************************
{
  // Allocate new memory entry to store the node Ids
  m_node = newEntry<Int>(m_nnodes, ValType::INT, VarType::SCALAR, NODEID_NAME);

  // Allocate new memory entry to store the coordinates of nodes
  m_coord = newEntry<Real>(m_nnodes, ValType::REAL, VarType::VECTOR,
                           COORD_NAME);

  // Allocate new memory entry to store the line element Ids
  m_linId = newEntry<Int>(m_nLins, ValType::INT, VarType::SCALAR, LINID_NAME);

  // Allocate new memory entry to store the triangle element Ids
  m_triId = newEntry<Int>(m_nTris, ValType::INT, VarType::SCALAR, TRIID_NAME);

  // Reserve capacity to store element connectivity and tags
  reserveElem();
}

void
GmshTxtMeshReader::addElem(vector<vector<Int>>& elemType, vector<Int>& nodes)
//******************************************************************************
//  Add new element (connectivity) to given element type container
//! \param[in]  elemType   Container type to add to (lines, triangles, etc)
//! \param[in]  nodes      Connectivity (node Ids) of the new element
//! \author J. Bakosi
//******************************************************************************
{
  try {
    elemType.push_back(nodes);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}

void
GmshTxtMeshReader::addElemTags(vector<vector<Int>>& elemType, vector<Int>& tags)
//******************************************************************************
//  Add new element tags to given element type
//! \param[in]  elemType  Container type to add to (lines, triangles, etc)
//! \param[in]  tags      Vector of tags to be added
//! \author J. Bakosi
//******************************************************************************
{
  try {
    elemType.push_back(tags);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}

void
GmshTxtMeshReader::reserveElem()
//******************************************************************************
//  Add new element
//! \param[in]  n  Desired new capacity to store n elements with their tags
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_lin.reserve(m_nLins);
    m_tri.reserve(m_nTris);
    m_linTag.reserve(m_nLins);
    m_triTag.reserve(m_nTris);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}
