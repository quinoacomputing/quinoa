//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 05:53:49 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtMeshWriter class declaration
  \details   GmshTxtMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshWriter_h
#define GmshTxtMeshWriter_h

#include <string>

using namespace std;

#include <MeshWriter.h>

namespace Quinoa {

//! GmshTxtMeshWriter : MeshWriter
class GmshTxtMeshWriter : private MeshWriter {

  public:
    //! Constructor
    GmshTxtMeshWriter(string filename, UnsMesh* mesh, Memory* memory) :
      MeshWriter(filename, mesh, memory) {}

    //! Destructor, default compiler generated
    ~GmshTxtMeshWriter() = default;

    //! Write Gmsh mesh to file
    void write();

  private:
    //! Don't permit copy constructor
    GmshTxtMeshWriter(const GmshTxtMeshWriter&) = delete;
    //! Don't permit copy assigment
    GmshTxtMeshWriter& operator=(const GmshTxtMeshWriter&) = delete;
    //! Don't permit move constructor
    GmshTxtMeshWriter(GmshTxtMeshWriter&&) = delete;
    //! Don't permit move assigment
    GmshTxtMeshWriter& operator=(GmshTxtMeshWriter&&) = delete;

    //! Write mandatory "$MeshFormat--$EndMeshFormat" section
    void writeMeshFormat();

    //! Write "$Nodes--$EndNodes" section
    void writeNodes();

    //! Write "$Elements--$EndElements" section
    void writeElements();

    //! Write "$PhysicalNames--$EndPhysicalNames" section
    void writePhysicalNames();
};

} // namespace Quinoa

#endif // GmshTxtMeshWriter_h
