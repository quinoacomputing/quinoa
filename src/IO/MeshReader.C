//******************************************************************************
/*!
  \file      src/Mesh/MeshReader.C
  \author    J. Bakosi
  \date      Tue 11 Sep 2012 06:25:39 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh reader class definition
  \details   Mesh reader class definition
*/
//******************************************************************************

#include <MeshReader.h>
#include <MeshException.h>

using namespace Quinoa;

MeshReader::MeshReader(string filename, UnsMesh* mesh, Memory* memory) :
              m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Open Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  // Open mesh file
  m_inMesh.open(m_filename, ifstream::in);
  if (!m_inMesh.good()) throw MeshException(FATAL, FAILED_OPEN, m_filename);
}

MeshReader::~MeshReader()
//******************************************************************************
//  Destructor: Close Gmsh mesh file
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.close();
  if (m_inMesh.fail()) throw MeshException(WARNING, FAILED_CLOSE, m_filename);
}
