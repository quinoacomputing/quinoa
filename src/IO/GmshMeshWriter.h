//******************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 08:43:57 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Gmsh mesh writer class declaration
  \details   Gmsh mesh writer class declaration. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
//******************************************************************************
#ifndef GmshMeshWriter_h
#define GmshMeshWriter_h

#include <string>

#include <Writer.h>
#include <UnsMesh.h>
#include <GmshMeshIO.h>

namespace quinoa {

//! \brief GmshMeshWriter : Writer
//! \details Mesh writer class facilitating writing a mesh to a file readable by
//!   the Gmsh mesh generator: http://geuz.org/gmsh.
class GmshMeshWriter : public tk::Writer {

  public:
    //! Constructor
    explicit GmshMeshWriter( const std::string& filename,
                             const UnsMesh& mesh,
                             GmshFileType type = GmshFileType::BINARY,
                             tk::real version = 2.2,
                             int datasize = sizeof(double) );

    //! Write Gmsh mesh to file
    void write() override;

  private:
    //! Write "$Nodes--$EndNodes" section
    void writeNodes();

    //! Write "$Elements--$EndElements" section
    void writeElements();

    //! Write "$PhysicalNames--$EndPhysicalNames" section
    void writePhysicalNames();

    //! \brief Mesh ASCII type query
    //! \return true if member variable m_type indicates an ASCII mesh format
    bool isASCII() const {
      return m_type == GmshFileType::ASCII ? true : false;
    }
    //! \brief Mesh binary type query
    //! \return true if member variable m_type indicates an binary mesh format
    bool isBinary() const {
      return m_type == GmshFileType::BINARY ? true : false;
    }

    //! Write element block: element ids, tags, and connectivity (node list)
    void writeElemBlock( GmshElemType type,
                         const std::vector< int >& id,
                         const std::vector< std::vector< int > >& tag,
                         const std::vector< std::vector< int > >& inpoel );

    const UnsMesh& m_mesh;              //!< Mesh object
    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary
};

} // quinoa::

#endif // GmshMeshWriter_h
