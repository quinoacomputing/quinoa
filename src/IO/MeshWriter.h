//******************************************************************************
/*!
  \file      src/IO/MeshWriter.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:11:15 AM KST
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

  public:
    //! Constructor: Acquire mesh file handle
    MeshWriter(string filename, UnsMesh* mesh, Memory* memory);

    //! Destructor: Release mesh file handle
    virtual ~MeshWriter();

    //! Interface for write mesh
    virtual void write() = 0;

  protected:
    //! Mesh file name
    string m_filename;

    //! Mesh file output stream
    ofstream m_outMesh;

    //! Mesh object pointer
    UnsMesh* m_mesh;

    //! Memory object pointer
    Memory* m_memory;

  private:
    //! Don't permit copy operator
    MeshWriter(const MeshWriter&);
    //! Don't permit assigment operator
    MeshWriter& operator=(const MeshWriter&);
};

} // namespace Quinoa

#endif // MeshWriter_h
