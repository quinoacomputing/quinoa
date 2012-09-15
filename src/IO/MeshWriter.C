//******************************************************************************
/*!
  \file      src/IO/MeshWriter.C
  \author    J. Bakosi
  \date      Sat 15 Sep 2012 02:13:47 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh writer class definition
  \details   Mesh writer class definition
*/
//******************************************************************************

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
    throw IOException(ExceptType::FATAL, IOExceptType::FAILED_OPEN, m_filename);
}

MeshWriter::~MeshWriter()
//******************************************************************************
//  Destructor: Release mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh.close();
  if (m_outMesh.fail())
    throw IOException(ExceptType::WARNING,
                      IOExceptType::FAILED_CLOSE,
                      m_filename);
}
