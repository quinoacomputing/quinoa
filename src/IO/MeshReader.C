//******************************************************************************
/*!
  \file      src/IO/MeshReader.C
  \author    J. Bakosi
  \date      Fri Sep 14 17:50:42 2012
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
    throw IOException(ExceptType::FATAL, IOExceptType::FAILED_OPEN, m_filename);
}

MeshReader::~MeshReader()
//******************************************************************************
//  Destructor: Release mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_inMesh.close();
  if (m_inMesh.fail())
    throw IOException(ExceptType::WARNING,
                      IOExceptType::FAILED_CLOSE,
                      m_filename);
}
