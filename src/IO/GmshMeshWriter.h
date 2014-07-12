//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \date      Wed Apr 23 11:14:58 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     GmshMeshWriter class declaration
  \details   GmshMeshWriter class declaration
*/
//******************************************************************************
#ifndef GmshMeshWriter_h
#define GmshMeshWriter_h

#include <string>

#include <Writer.h>
#include <UnsMesh.h>
#include <GmshMeshIO.h>

namespace quinoa {

//! GmshMeshWriter : Writer
class GmshMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit GmshMeshWriter( const std::string& filename,
                             UnsMesh& mesh,
                             GmshFileType type = GmshFileType::BINARY,
                             tk::real version = 2.2,
                             int datasize = sizeof(double) );

    //! Destructor, default compiler generated
    ~GmshMeshWriter() noexcept override = default;

    //! Write Gmsh mesh to file
    void write() override;

  private:
    //! Don't permit copy constructor
    GmshMeshWriter(const GmshMeshWriter&) = delete;
    //! Don't permit copy assigment
    GmshMeshWriter& operator=(const GmshMeshWriter&) = delete;
    //! Don't permit move constructor
    GmshMeshWriter(GmshMeshWriter&&) = delete;
    //! Don't permit move assigment
    GmshMeshWriter& operator=(GmshMeshWriter&&) = delete;

    //! Write "$Nodes--$EndNodes" section
    void writeNodes();

    //! Write "$Elements--$EndElements" section
    void writeElements();

    //! Write "$PhysicalNames--$EndPhysicalNames" section
    void writePhysicalNames();

    // Get mesh type
    bool isASCII() const {
      return m_type == GmshFileType::ASCII ? true : false;
    }
    bool isBinary() const {
      return m_type == GmshFileType::BINARY ? true : false;
    }

    //! Write element block: element ids, tags, and connectivity (node list)
    void writeElemBlock( GmshElemType type, std::vector< int >& id,
                         std::vector< std::vector< int > >& tag,
                         std::vector< std::vector< int > >& inpoel );

    UnsMesh& m_mesh;                    //!< Mesh object
    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary
};

} // quinoa::

#endif // GmshMeshWriter_h
