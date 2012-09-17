//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 08:54:23 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <GmshMeshWriter.h>
#include <IOException.h>

using namespace Quinoa;

void
GmshMeshWriter::write()
//******************************************************************************
//  Write Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  // Write out mandatory "$MeshFormat" section
  writeMeshFormat();

//   // Keep reading in sections until end of file
//   while (!m_inMesh.eof()) {
//     string s;
//     getline(m_inMesh, s);
//     if (s=="$Nodes") readNodes();
//     else if (s=="$Elements") readElements();
//     else if (s=="$PhysicalNames") readPhysicalNames();
//   }
// 
//   // Clear failbit triggered by eof, so close() won't throw a false FAILED_CLOSE
//   m_inMesh.clear();
}

void
GmshMeshWriter::writeMeshFormat()
//******************************************************************************
//  Write mandatory "$MeshFormat" section
//! \author J. Bakosi
//******************************************************************************
{
  string s;

  // Write beginning of header: $MeshFormat
  m_outMesh << "$MeshFormat\n";
  if (m_outMesh.bad())
    throw IOException(ExceptType::FATAL,
                      IOExceptType::FAILED_WRITE, 
                      m_filename);

  // Write "version-number file-type data-size"
  m_outMesh << m_mesh->getVersion() << " "
            << m_mesh->getType() << " "
            << m_mesh->getDatasize() << "\n";
  if (m_outMesh.bad())
    throw IOException(ExceptType::FATAL,
                      IOExceptType::FAILED_WRITE, 
                      m_filename);

  // Write end of header: $EndMeshFormat
  m_outMesh << "$EndMeshFormat" << endl;
  if (m_outMesh.bad())
    throw IOException(ExceptType::FATAL,
                      IOExceptType::FAILED_WRITE, 
                      m_filename);
}
