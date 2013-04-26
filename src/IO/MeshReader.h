//******************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:51:10 2013
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
    explicit MeshReader(const string filename,
                        UnsMesh* const mesh,
                        Memory* const memory);

    //! Destructor: Release mesh file handle
    virtual ~MeshReader() noexcept;

    //! Interface for read mesh
    virtual void read() = 0;

    const string m_filename;            //!< Mesh file name
    UnsMesh* const m_mesh;              //!< Mesh object pointer
    Memory* const m_memory;             //!< Memory object pointer

    //! Mesh file input stream
    ifstream m_inMesh;

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
