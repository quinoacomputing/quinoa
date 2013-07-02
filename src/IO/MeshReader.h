//******************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \author    J. Bakosi
  \date      Tue Jul  2 15:25:37 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshReader base class declaration
  \details   MeshReader base class declaration
*/
//******************************************************************************
#ifndef MeshReader_h
#define MeshReader_h

#include <fstream>

#include <UnsMesh.h>
#include <Memory.h>

namespace Quinoa {

//! MeshReader base class
class MeshReader {

  protected:
    //! Constructor: Acquire mesh file handle
    explicit MeshReader(const std::string filename,
                        UnsMesh* const mesh,
                        Memory* const memory);

    //! Destructor: Release mesh file handle
    virtual ~MeshReader() noexcept;

    //! Interface for read mesh
    virtual void read() = 0;

    const std::string m_filename;            //!< Mesh file name
    UnsMesh* const m_mesh;                   //!< Mesh object pointer
    Memory* const m_memory;                  //!< Memory object pointer

    //! Mesh file input stream
    std::ifstream m_inMesh;

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
