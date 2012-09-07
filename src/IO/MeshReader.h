//******************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 04:26:00 PM MDT
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

#include <Mesh.h>
#include <Exception.h>

namespace Quinoa {

//! MeshReader base class
class MeshReader {

  public:
    //! Constructor
    MeshReader(string filename) : m_filename(filename) {}

    //! Destructor
    ~MeshReader() {};

    //! Interface for open
    virtual void open() = 0;

    //! Interface for close
    virtual void close() = 0;

    //! Interface for read mesh
    virtual void read(UnsMesh* mesh) = 0;

  protected:
    //! Mesh file name
    string m_filename;
    //! Mesh file input stream
    ifstream m_mesh;

  private:
    //! Don't permit copy operator
    MeshReader(const MeshReader&);
    //! Don't permit assigment operator
    MeshReader& operator=(const MeshReader&);
};

} // namespace Quinoa

#endif // MeshReader_h
