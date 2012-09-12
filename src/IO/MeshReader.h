//******************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 03:54:13 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshReader base class declaration
  \details   MeshReader base class declaration
*/
//******************************************************************************
#ifndef MeshReader_h
#define MeshReader_h

#include <string>
#include <fstream>

using namespace std;

#include <UnsMesh.h>
#include <Exception.h>
#include <Memory.h>

namespace Quinoa {

//! MeshReader base class
class MeshReader {

  public:
    //! Constructor: Acquire mesh file handle
    MeshReader(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release mesh file handle
    virtual ~MeshReader();

    //! Interface for read mesh
    virtual void read() = 0;

  protected:
    //! Mesh file name
    string m_filename;

    //! Mesh file input stream
    ifstream m_inMesh;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy operator
    MeshReader(const MeshReader&);
    //! Don't permit assigment operator
    MeshReader& operator=(const MeshReader&);
};

} // namespace Quinoa

#endif // MeshReader_h
