//******************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 08:25:36 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshReader base class declaration
  \details   MeshReader base class declaration
*/
//******************************************************************************
#ifndef MeshReader_h
#define MeshReader_h

#include <fstream>

using namespace std;

#include <UnsMesh.h>
#include <Memory.h>

namespace Quinoa {

//! MeshReader base class
class MeshReader {

  protected:
    //! Constructor: Acquire mesh file handle
    MeshReader(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release mesh file handle
    ~MeshReader();

    //! Interface for read mesh
    virtual void read() = 0;

    //! Mesh file name
    string m_filename;

    //! Mesh file input stream
    ifstream m_inMesh;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy constructor
    MeshReader(const MeshReader&) = delete;
    //! Don't permit copy assigment
    MeshReader& operator=(const MeshReader&) = delete;
    //! Don't permit move constructor
    MeshReader(MeshReader&&) = delete;
    //! Don't permit move assigment
    MeshReader& operator=(MeshReader&&) = delete;
};

} // namespace Quinoa

#endif // MeshReader_h
