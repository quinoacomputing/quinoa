//******************************************************************************
/*!
  \file      src/IO/MeshReader.C
  \author    J. Bakosi
  \date      Mon May  6 17:24:37 2013
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
  ErrChk(m_inMesh.good(), FATAL, "Failed to open file: " + m_filename);
}

MeshReader::~MeshReader() noexcept
//******************************************************************************
//  Destructor: Release mesh file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_inMesh.close();
    ErrChk(!m_inMesh.fail(), WARNING, "Failed to close file: " + m_filename);

  } // emit only a warning on error
    catch (Exception& e) {
      e.echo("WARNING");
    }
    catch (exception& e) {
      cout << ">>> std::exception in MeshReader destructor: " << e.what()
           << endl;
    }
    catch (...) {
      cout << ">>> UNKNOWN EXCEPTION in MeshReader destructor" << endl;
    }
}
