//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 05:42:22 AM KST
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

    //! Destructor
    ~GmshMeshWriter() {};

    //! Write Gmsh mesh to file
    void write();

  private:
    //! Don't permit copy operator
    GmshMeshWriter(const GmshMeshWriter&);
    //! Don't permit assigment operator
    GmshMeshWriter& operator=(const GmshMeshWriter&);
};

} // namespace Quinoa

#endif // GmshMeshWriter_h
