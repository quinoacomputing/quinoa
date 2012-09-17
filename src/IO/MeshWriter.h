//******************************************************************************
/*!
  \file      src/IO/MeshWriter.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 05:57:47 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshWriter base class declaration
  \details   MeshWriter base class declaration
*/
//******************************************************************************
#ifndef MeshWriter_h
#define MeshWriter_h

#include <string>
#include <fstream>

using namespace std;

#include <UnsMesh.h>
#include <Memory.h>

namespace Quinoa {

//! MeshWriter base class
class MeshWriter {

  protected:
    //! Constructor: Acquire mesh file handle
    MeshWriter(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release mesh file handle
    virtual ~MeshWriter();

    //! Interface for write mesh
    virtual void write() = 0;

    //! Mesh file name
    string m_filename;

    //! Mesh file output stream
    ofstream m_outMesh;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy constructor
    MeshWriter(const MeshWriter&) = delete;
    //! Don't permit copy assigment
    MeshWriter& operator=(const MeshWriter&) = delete;
    //! Don't permit move constructor
    MeshWriter(MeshWriter&&) = delete;
    //! Don't permit move assigment
    MeshWriter& operator=(MeshWriter&&) = delete;
};

} // namespace Quinoa

#endif // MeshWriter_h
