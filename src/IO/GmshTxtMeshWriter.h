//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.h
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 02:19:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtMeshWriter class declaration
  \details   GmshTxtMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshWriter_h
#define GmshTxtMeshWriter_h

#include <string>

#include <Writer.h>
#include <GmshMesh.h>

namespace quinoa {

//! GmshTxtMeshWriter : Writer
class GmshTxtMeshWriter : public Writer {

  public:
    //! Constructor
    explicit GmshTxtMeshWriter(const std::string& filename,
                              GmshMesh* const mesh) :
      Writer(filename),
      m_mesh(mesh) {}

    //! Destructor, default compiler generated
    ~GmshTxtMeshWriter() noexcept override = default;

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

    GmshMesh* const m_mesh;         //!< Mesh object pointer
};

} // namespace quinoa

#endif // GmshTxtMeshWriter_h
