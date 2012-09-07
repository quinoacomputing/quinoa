//******************************************************************************
/*!
  \file      src/IO/GmshReader.h
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 03:35:39 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshReader_h
#define GmshReader_h

#include <UnsMesh.h>
#include <MeshReader.h>
#include <GmshException.h>

namespace Quinoa {

//! GmshReader : MeshReader
class GmshReader : MeshReader {

  public:
    //! Constructor
    GmshReader(string filename) : MeshReader(filename) {}

    //! Destructor
    ~GmshReader() {};

    //! Interface for open
    virtual void open();

    //! Interface for close
    virtual void close();

    //! Interface for read
    virtual void read(Mesh* mesh);

  private:
    //! Don't permit copy operator
    GmshReader(const GmshReader&);
    //! Don't permit assigment operator
    GmshReader& operator=(const GmshReader&);
};

} // namespace Quinoa

#endif // GmshReader_h
