//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.C
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:41:09 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh mesh writer class definition
  \details   Gmsh mesh writer class definition
*/
//******************************************************************************

#include <GmshMeshWriter.h>

using namespace Quinoa;

void
GmshMeshWriter::write()
//******************************************************************************
//  Write Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
//   // Read in mandatory "$MeshFormat" section
//   readMeshFormat();
// 
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

