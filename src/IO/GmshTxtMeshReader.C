//******************************************************************************
/*!
  \file      src/Mesh/GmshTxtMeshReader.C
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 07:23:39 AM MDT
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
GmshTxtMeshReader::read()
//******************************************************************************
//  Public interface for read Gmsh mesh
//! \author J. Bakosi
//******************************************************************************
{
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
}

void
GmshTxtMeshReader::echoElemSets()
//******************************************************************************
//  Echo element tags and connectivity in all element sets
//! \author J. Bakosi
//******************************************************************************
{
  // Return and throw warning if there are no element sets
  if (!m_elemsets)
    throw MeshException(ExceptType::WARNING, MeshExceptType::EMPTY_SET);

  // Echo all element sets
  for (Int k=1; k<=m_elemsets; k++) {
    cout << "* Element set: " << k << endl << endl;
    stringstream ess;
    ess << k;
    Int* elmtype = m_memory->getPtr<Int>(ELEMTYPE_NAME+ess.str());
    pair<size_t,Int*> elmEntry = m_memory->getNumPtr<Int>(ELEMID_NAME+ess.str());
    size_t num = elmEntry.first;
    Int* element = elmEntry.second;

    // elm-number elm-type number-of-tags < tag > ... node-number-list
    for (size_t i=0; i<num; i++) {
      cout << "  " << element[i] << " " << elmtype[i] << " {";
      copy(m_tag[i].begin(),m_tag[i].end()-1,ostream_iterator<Int>(cout,", "));
      cout << m_tag[i].back()-1 << "} {";
      copy(m_elem[i].begin(),m_elem[i].end()-1,ostream_iterator<Int>(cout,", "));
      cout << m_elem[i].back()-1 << "}" << endl;
    }
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
GmshTxtMeshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  Int num;
  m_inMesh >> num;

  // Increase number of node sets
  stringstream nss;
  nss << ++m_nodesets;

  // Add new MeshSet entry for the node ids
  Int* node = newEntry<Int>(num,
                            ValType::INT,
                            VarType::SCALAR,
                            NODEID_NAME+nss.str());
  // Add new MeshSet entry for the coordinates of the node
  Real* coord = newEntry<Real>(num,
                               ValType::REAL,
                               VarType::VECTOR,
                               COORD_NAME+nss.str());

  // Read in node ids and coordinates
  for (Int i=0; i<num; i++) {
    // node-number x-coord y-coord z-coord
    m_inMesh >> node[i] >> coord[3*i] >> coord[3*i+1] >> coord[3*i+2];
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
GmshTxtMeshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  Int num;
  m_inMesh >> num;

  // Increase number of element sets
  stringstream ess;
  ess << ++m_elemsets;

  // Add new MeshSet entry for the element ids
  Int* elmId = newEntry<Int>(num,
                             ValType::INT,
                             VarType::SCALAR,
                             ELEMID_NAME+ess.str());
  // Add new MeshSet entry for the element types
  Int* elmType = newEntry<Int>(num,
                               ValType::INT,
                               VarType::SCALAR,
                               ELEMTYPE_NAME+ess.str());

  // Reserve capacity to store connectivity and tags
  reserveElem(num);

  // Read in element ids, tags, and element connectivity (node list)
  for (Int i=0; i<num; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    Int ntags;
    m_inMesh >> elmId[i] >> elmType[i] >> ntags;

    // Read and add tags
    vector<Int> tags(ntags,0);
    for (Int j=0; j<ntags; j++) m_inMesh >> tags[j];
    addElemTags(tags);

    // Read and add element node list
    auto it = GmshElemNodes.find(elmType[i]);
    if (it==GmshElemNodes.end())
      throw MeshException(ExceptType::FATAL,
            MeshExceptType::BAD_ELEMENT,
            m_filename);
    Int nnodes = it->second;
    vector<Int> nodes(nnodes,0);
    for (Int j=0; j<nnodes; j++) m_inMesh >> nodes[j];
    addElem(nodes);
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
GmshTxtMeshReader::addElem(vector<int>& nodes)
//******************************************************************************
//  Add new element
//! \param[in]  nodes  Vector of node ids (i.e. connectivity) of the new element
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_elem.push_back(nodes);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}

void
GmshTxtMeshReader::addElemTags(vector<Int>& tags)
//******************************************************************************
//  Add new element tags
//! \param[in]  tags  Vector of tags to be added
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_tag.push_back(tags);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}

void
GmshTxtMeshReader::reserveElem(vector<vector<Int>>::size_type n)
//******************************************************************************
//  Add new element
//! \param[in]  n  Desired new capacity to store n elements with their tags
//! \author J. Bakosi
//******************************************************************************
{
  try {
    m_elem.reserve(n);
    m_tag.reserve(n);
  } catch (bad_alloc& ba) {
    throw MemoryException(ExceptType::FATAL, MemExceptType::BAD_ALLOC);
  }
}

