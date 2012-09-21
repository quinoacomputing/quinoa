//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshReader2D.h
  \author    J. Bakosi
  \date      Fri 21 Sep 2012 09:37:54 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshReader2D_h
#define GmshTxtMeshReader2D_h

#include <GmshTxtMeshReader.h>

namespace Quinoa {

//! GmshTxtMeshReader2D : GmshTxtMeshReader
class GmshTxtMeshReader2D : GmshTxtMeshReader {

  public:
    //! Constructor
    GmshTxtMeshReader2D(string filename, UnsMesh* mesh, Memory* memory) :
      GmshTxtMeshReader(filename, mesh, memory) {}

    //! Destructor
    ~GmshTxtMeshReader2D() = default;

    //! Public interface for read 2D txt Gmsh mesh
    void read();

  private:
    //! Don't permit copy constructor
    GmshTxtMeshReader2D(const GmshTxtMeshReader2D&) = delete;
    //! Don't permit copy assigment
    GmshTxtMeshReader2D& operator=(const GmshTxtMeshReader2D&) = delete;
    //! Don't permit move constructor
    GmshTxtMeshReader2D(GmshTxtMeshReader2D&&) = delete;
    //! Don't permit move assigment
    GmshTxtMeshReader2D& operator=(GmshTxtMeshReader2D&&) = delete;

    //! Compress mesh data from STL to memory entries
    void compress();
};

} // namespace Quinoa

#endif // GmshTxtMeshReader2D_h
