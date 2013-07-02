//******************************************************************************
/*!
  \file      src/IO/MeshWriter.C
  \author    J. Bakosi
  \date      Tue Jul  2 15:28:13 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mesh writer class definition
  \details   Mesh writer class definition
*/
//******************************************************************************

#include <iostream>

#include <MeshWriter.h>

using namespace Quinoa;

MeshWriter::MeshWriter(std::string filename, UnsMesh* mesh, Memory* memory) :
              m_filename(filename), m_mesh(mesh), m_memory(memory)
//******************************************************************************
//  Constructor: Acquire mesh file handle
//! \author J. Bakosi
//******************************************************************************
{
  m_outMesh.open(m_filename, std::ofstream::out);
  ErrChk(m_outMesh.good(), ExceptType::FATAL,
         "Failed to open file: " + m_filename);
}

MeshWriter::~MeshWriter() noexcept
//******************************************************************************
//  Destructor: Release mesh file handle
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_outMesh.close();
    ErrChk(!m_outMesh.fail(), ExceptType::WARNING,
           "Failed to close file: " + m_filename);

  } // emit only a warning on error
    catch (Exception& e) {
      e.echo("WARNING");
    }
    catch (std::exception& e) {
      std::cout << ">>> std::exception in MeshWriter destructor: " << e.what()
                << std::endl;
    }
    catch (...) {
      std::cout << ">>> UNKNOWN EXCEPTION in MeshWriter destructor"
                << std::endl;
    }
}
