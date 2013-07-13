//******************************************************************************
/*!
  \file      src/IO/GmshTxtMeshWriter.h
  \author    J. Bakosi
  \date      Fri 12 Jul 2013 09:56:18 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshTxtMeshWriter class declaration
  \details   GmshTxtMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshTxtMeshWriter_h
#define GmshTxtMeshWriter_h

#include <string>

#include <Writer.h>
#include <UnsMesh.h>

namespace Quinoa {

//! GmshTxtMeshWriter : Writer
class GmshTxtMeshWriter : public Writer {

  public:
    //! Constructor
    explicit GmshTxtMeshWriter(const std::string filename,
                              UnsMesh* const mesh) :
      Writer(filename),
      m_mesh(mesh) {}

    //! Destructor, default compiler generated
    virtual ~GmshTxtMeshWriter() noexcept = default;

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

    UnsMesh* const m_mesh;         //!< Mesh object pointer
};

} // namespace Quinoa

#endif // GmshTxtMeshWriter_h
