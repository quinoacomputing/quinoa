//******************************************************************************
/*!
  \file      src/IO/MeshReader.C
  \author    J. Bakosi
  \date      Wed 01 May 2013 09:30:36 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh reader class definition
  \details   Mesh reader class definition
*/
//******************************************************************************

#include <iostream>

#include <MeshReader.h>

using namespace Quinoa;

MeshReader::MeshReader(const string filename,
                       UnsMesh* const mesh,
                       Memory* const memory) :
  m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Acquire mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.open(m_filename, ifstream::in);
  if (!m_inMesh.good())
    throw Exception( FATAL, "Failed to open file: " + m_filename);
}

MeshReader::~MeshReader() noexcept
//******************************************************************************
//  Destructor: Release mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.close();

  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
  if (m_inMesh.fail())
    cout << "WARNING: Failed to close file: " << m_filename << endl;
}
