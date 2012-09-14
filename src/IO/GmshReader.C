//******************************************************************************
/*!
  \file      src/Mesh/GmshReader.C
  \author    J. Bakosi
  \date      Fri Sep 14 16:20:07 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh reader class definition
  \details   Gmsh mesh reader class definition
*/
//******************************************************************************

#include <fstream>
#include <sstream>
#include <iostream>

#include <QuinoaTypes.h>
#include <GmshReader.h>
#include <MeshException.h>
#include <Memory.h>

using namespace Quinoa;

void
GmshReader::read()
//******************************************************************************
//  Read Gmsh mesh from file
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
GmshReader::readMeshFormat()
//******************************************************************************
//  Read mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  string s;

  // Read in beginning of header: $MeshFormat
  getline(m_inMesh, s);
  if (s!="$MeshFormat") throw MeshException(FATAL, BAD_FORMAT, m_filename);

  // Read in "version-number file-type data-size"
  Real version;
  Int type, datasize;
  m_inMesh >> version >> type >> datasize;
  if ((version!=2.2 && version!=2.0) || type!=0 || datasize!=sizeof(Real))
    throw MeshException(FATAL, BAD_FORMAT, m_filename);
  getline(m_inMesh, s);  // finish reading the line

  // Read in end of header: $EndMeshFormat
  getline(m_inMesh, s);
  if (s!="$EndMeshFormat") throw MeshException(FATAL, BAD_FORMAT, m_filename);
}

void
GmshReader::readNodes()
//******************************************************************************
//  Read "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of nodes in this node set
  Int num;
  m_inMesh >> num;

  // Increase number of node sets
  Int nodesets = m_mesh->addNodeSet();
  stringstream nss;
  nss << nodesets;

  // Add new MeshSet entry for the node ids
  Int* node = m_mesh->newEntry<Int>(num,
                                    ValType::INT,
                                    VarType::SCALAR,
                                    NODEID_NAME+nss.str());
  // Add new MeshSet entry for the coordinates of the node
  Real* coord = m_mesh->newEntry<Real>(num,
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
  if (s!="$EndNodes") throw MeshException(FATAL, BAD_FORMAT, m_filename);
}

void
GmshReader::readElements()
//******************************************************************************
//  Read "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  // Read in number of elements in this element set
  Int num;
  m_inMesh >> num;

  // Increase number of element sets
  Int elemsets = m_mesh->addElemSet();
  stringstream ess;
  ess << elemsets;

  // Add new MeshSet entry for the element ids
  Int* element = m_mesh->newEntry<Int>(num,
                                       ValType::INT,
                                       VarType::SCALAR,
                                       ELEMID_NAME+ess.str());
  // Add new MeshSet entry for the element types
  Int* elmtype = m_mesh->newEntry<Int>(num,
                                       ValType::INT,
                                       VarType::SCALAR,
                                       ELEMTYPE_NAME+ess.str());

  // Reserve capacity to store connectivity and tags
  m_mesh->reserveElem(num);

  // Read in element ids, tags, and node list
  for (Int i=0; i<num; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    Int ntags;
    m_inMesh >> element[i] >> elmtype[i] >> ntags;

    // Read and add tags
    vector<Int> tags(ntags,0);
    for (Int j=0; j<ntags; j++) {
      m_inMesh >> tags[j];
    }
    m_mesh->addElemTags(tags);

    // Read and add element node list
    auto it = GmshElemNodes.find(elmtype[i]);
    if (it==GmshElemNodes.end()) {
      throw MeshException(FATAL, BAD_ELEMENT, m_filename);
    }
    Int nnodes = it->second;
    vector<Int> nodes(nnodes,0);
    for (Int j=0; j<nnodes; j++) {
      m_inMesh >> nodes[j];
    }
    m_mesh->addElem(nodes);
  }
  string s;
  getline(m_inMesh, s);  // finish reading the last line

  // Read in end of header: $EndNodes
  getline(m_inMesh, s);
  if (s!="$EndElements") throw MeshException(FATAL, BAD_FORMAT, m_filename);
}

void
GmshReader::readPhysicalNames()
//******************************************************************************
//  Read "$PhysicalNames--$EndPhysicalNames" section
//! \author J. Bakosi
//******************************************************************************
{
}
