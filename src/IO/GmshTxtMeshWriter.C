//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.C
  \author    J. Bakosi
  \date      Tue Jul  2 15:30:45 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <iterator>
#include <iomanip>

#include <GmshTxtMeshWriter.h>

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
  // Write beginning of header: $MeshFormat
  m_outMesh << "$MeshFormat\n";
  ErrChk(!m_outMesh.bad(), ExceptType::FATAL,
         "Failed to write to file: " + m_filename);

  // Write "version-number file-type data-size"
  m_outMesh << m_mesh->getVersion() << " "
            << m_mesh->getType() << " "
            << m_mesh->getDatasize() << "\n";
  ErrChk(!m_outMesh.bad(), ExceptType::FATAL,
         "Failed to write to file: " + m_filename);

  // Write end of header: $EndMeshFormat
  m_outMesh << "$EndMeshFormat" << std::endl;
  ErrChk(!m_outMesh.bad(), ExceptType::FATAL,
         "Failed to write to file: " + m_filename);
}

void
GmshTxtMeshWriter::writeNodes()
//******************************************************************************
//  Write "$Nodes--$EndNodes" section
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh << "$Nodes" << std::endl;

  // Get number of nodes, and pointers node ids and coordinates
  int nnodes = m_mesh->getNnodes();
  int* nodeId = m_mesh->getNodeId();
  real* coord = m_mesh->getCoord();

  // Write out number of nodes
  m_outMesh << nnodes << std::endl;

  // Write node ids and coordinates
  for (int i=0; i<nnodes; ++i) {
    // node-number x-coord y-coord z-coord
    int i3 = i*3;
    m_outMesh << nodeId[i] << " " << std::setprecision(16)
              << coord[i3] << " " << coord[i3+1] << " " << coord[i3+2]
              << std::endl;
  }

  m_outMesh << "$EndNodes" << std::endl;
}

void
GmshTxtMeshWriter::writeElements()
//******************************************************************************
//  Write "$Elements--$EndElements" section
//! \author J. Bakosi
//******************************************************************************
{
  using ST = std::vector<std::vector<int>>::size_type;

  m_outMesh << "$Elements" << std::endl;

  // Get pointers to the element ids, connectivities, and tags
  int* linId = m_mesh->getLineId();
  int* triId = m_mesh->getTriangleId();
  std::vector<std::vector<int>> linpoel = m_mesh->getLinpoel();
  std::vector<std::vector<int>> lintag = m_mesh->getLintag();
  std::vector<std::vector<int>> tinpoel = m_mesh->getTinpoel();
  std::vector<std::vector<int>> tritag = m_mesh->getTritag();

  // Get number of line and triangle elements
  ST nlines = linpoel.size();
  ST ntriangles = tinpoel.size();

  // Write out number of nodes
  m_outMesh << nlines+ntriangles << std::endl;

  // Write out line element ids, tags, and connectivity (node list)
  for (ST i=0; i<nlines; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outMesh << linId[i] << " " << 1 << " " << lintag[i].size() << " ";

    copy(lintag[i].begin(), lintag[i].end()-1,
         std::ostream_iterator<int>(m_outMesh," "));
    m_outMesh << lintag[i].back() << " ";

    copy(linpoel[i].begin(), linpoel[i].end()-1,
         std::ostream_iterator<int>(m_outMesh," "));
    m_outMesh << linpoel[i].back() << std::endl;
  }

  // Write out triangle element ids, tags, and connectivity (node list)
  for (ST i=0; i<ntriangles; i++) {
    // elm-number elm-type number-of-tags < tag > ... node-number-list
    m_outMesh << triId[i] << " " << 2 << " " << tritag[i].size() << " ";

    copy(tritag[i].begin(), tritag[i].end()-1,
         std::ostream_iterator<int>(m_outMesh," "));
    m_outMesh << tritag[i].back() << " ";

    copy(tinpoel[i].begin(), tinpoel[i].end()-1,
         std::ostream_iterator<int>(m_outMesh," "));
    m_outMesh << tinpoel[i].back() << std::endl;
  }

  m_outMesh << "$EndElements" << std::endl;
}
