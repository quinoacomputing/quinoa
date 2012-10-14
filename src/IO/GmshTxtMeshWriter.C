//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 06:34:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <iterator>
#include <iomanip>

#include <GmshTxtMeshWriter.h>
#include <IOException.h>

using namespace Quinoa;

void
GmshTxtMeshWriter::write()
//******************************************************************************
//  Write Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  // Write out mandatory "$MeshFormat" section
  writeMeshFormat();

  // Write sections
  writeNodes();
  writeElements();
  //writePhysicalNames();
}

void
GmshTxtMeshWriter::writeMeshFormat()
//******************************************************************************
//  Write mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  string s;

  // Write beginning of header: $MeshFormat
  m_outMesh << "$MeshFormat\n";
  if (m_outMesh.bad())
    throw IOException(FATAL, FAILED_WRITE, m_filename);

  // Write "version-number file-type data-size"
  m_outMesh << m_mesh->getVersion() << " "
            << m_mesh->getType() << " "
            << m_mesh->getDatasize() << "\n";
  if (m_outMesh.bad())
    throw IOException(FATAL, FAILED_WRITE,  m_filename);

  // Write end of header: $EndMeshFormat
  m_outMesh << "$EndMeshFormat" << endl;
  if (m_outMesh.bad())
    throw IOException(FATAL, FAILED_WRITE, m_filename);
}

void
GmshTxtMeshWriter::writeNodes()
//******************************************************************************
//  Write "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh << "$Nodes" << endl;

  // Get number of nodes, and pointers node ids and coordinates
  Int nnodes = m_mesh->getNnodes();
  Int* nodeId = m_mesh->getNodeId();
  Real* coord = m_mesh->getCoord();

  // Write out number of nodes
  m_outMesh << nnodes << endl;

  // Write node ids and coordinates
  for (Int i=0; i<nnodes; ++i) {
    // node-number x-coord y-coord z-coord
    Int i3 = i*3;
    m_outMesh << nodeId[i] << " " << setprecision(16)
              << coord[i3] << " " << coord[i3+1] << " " << coord[i3+2] << endl;
  }

  m_outMesh << "$EndNodes" << endl;
}

void
GmshTxtMeshWriter::writeElements()
//******************************************************************************
//  Write "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  typedef vector<vector<Int>>::size_type ST;

  m_outMesh << "$Elements" << endl;

  // Get pointers to the element ids, connectivities, and tags
  Int* linId = m_mesh->getLineId();
  Int* triId = m_mesh->getTriangleId();
  vector<vector<Int>> linpoel = m_mesh->getLinpoel();
  vector<vector<Int>> lintag = m_mesh->getLintag();
  vector<vector<Int>> tinpoel = m_mesh->getTinpoel();
  vector<vector<Int>> tritag = m_mesh->getTritag();

  // Get number of line and triangle elements
  ST nlines = linpoel.size();
  ST ntriangles = tinpoel.size();

  // Write out number of nodes
  m_outMesh << nlines+ntriangles << endl;

  // Write out line element ids, tags, and connectivity (node list)
  for (ST i=0; i<nlines; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outMesh << linId[i] << " " << 1 << " " << lintag[i].size() << " ";

    copy(lintag[i].begin(), lintag[i].end()-1,
         ostream_iterator<Int>(m_outMesh," "));
    m_outMesh << lintag[i].back() << " ";

    copy(linpoel[i].begin(), linpoel[i].end()-1,
         ostream_iterator<Int>(m_outMesh," "));
    m_outMesh << linpoel[i].back() << endl;
  }

  // Write out triangle element ids, tags, and connectivity (node list)
  for (ST i=0; i<ntriangles; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outMesh << triId[i] << " " << 2 << " " << tritag[i].size() << " ";

    copy(tritag[i].begin(), tritag[i].end()-1,
         ostream_iterator<Int>(m_outMesh," "));
    m_outMesh << tritag[i].back() << " ";

    copy(tinpoel[i].begin(), tinpoel[i].end()-1,
         ostream_iterator<Int>(m_outMesh," "));
    m_outMesh << tinpoel[i].back() << endl;
  }

  m_outMesh << "$EndElements" << endl;
}
