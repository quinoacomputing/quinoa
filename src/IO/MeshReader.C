//******************************************************************************
/*!
  \file      src/IO/MeshReader.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 10:22:42 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh reader class definition
  \details   Mesh reader class definition
*/
//******************************************************************************

#include <MeshReader.h>
#include <IOException.h>

using namespace Quinoa;

MeshReader::MeshReader(string filename, UnsMesh* mesh, Memory* memory) :
              m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Acquire mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.open(m_filename, ifstream::in);
  if (!m_inMesh.good())
    throw IOException(FATAL, FAILED_OPEN, m_filename);
}

MeshReader::~MeshReader()
//******************************************************************************
//  Destructor: Release mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.close();

  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_inMesh.fail())
    cerr << "WARNING: Failed to close file: " << m_filename << endl;
}
