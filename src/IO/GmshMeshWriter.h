//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 06:26:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshMeshWriter class declaration
  \details   GmshMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshMeshWriter_h
#define GmshMeshWriter_h

#include <string>

using namespace std;

#include <MeshWriter.h>

namespace Quinoa {

//! GmshMeshWriter : MeshWriter
class GmshMeshWriter : MeshWriter {

  public:
    //! Constructor
    GmshMeshWriter(string filename, UnsMesh* mesh, Memory* memory) :
      MeshWriter(filename, mesh, memory) {}

    //! Destructor, default compiler generated
    ~GmshMeshWriter() = default;

    //! Write Gmsh mesh to file
    void write();

  private:
    //! Don't permit copy constructor
    GmshMeshWriter(const GmshMeshWriter&) = delete;
    //! Don't permit copy assigment
    GmshMeshWriter& operator=(const GmshMeshWriter&) = delete;
    //! Don't permit move constructor
    GmshMeshWriter(GmshMeshWriter&&) = delete;
    //! Don't permit move assigment
    GmshMeshWriter& operator=(GmshMeshWriter&&) = delete;
};

} // namespace Quinoa

#endif // GmshMeshWriter_h
