//******************************************************************************
/*!
  \file      src/IO/MeshWriter.C
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 10:05:49 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh writer class definition
  \details   Mesh writer class definition
*/
//******************************************************************************

#include <iostream>

#include <MeshWriter.h>
#include <IOException.h>

using namespace Quinoa;

MeshWriter::MeshWriter(string filename, UnsMesh* mesh, Memory* memory) :
              m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Acquire mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh.open(m_filename, ofstream::out);
  if (!m_outMesh.good())
    throw IOException(FATAL, IO_FAILED_OPEN, m_filename);
}

MeshWriter::~MeshWriter()
//******************************************************************************
//  Destructor: Release mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh.close();
  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_outMesh.fail())
    cerr << "WARNING: Failed to close file: " << m_filename << endl;
}
