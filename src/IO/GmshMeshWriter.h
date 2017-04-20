// *****************************************************************************
/*!
  \file      src/IO/GmshMeshWriter.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Gmsh mesh writer class declaration
  \details   Gmsh mesh writer class declaration. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
// *****************************************************************************
#ifndef GmshMeshWriter_h
#define GmshMeshWriter_h

#include <iosfwd>
#include <cstddef>
#include <vector>

#include "Types.h"
#include "Writer.h"
#include "GmshMeshIO.h"

namespace tk {

class UnsMesh;

//! Gmsh mesh writer
//! \details Mesh writer class facilitating writing a mesh to a file readable by
//!   the Gmsh mesh generator: http://geuz.org/gmsh.
class GmshMeshWriter : public Writer {

  public:
    //! Constructor
    explicit GmshMeshWriter( const std::string& filename,
                             GmshFileType type = GmshFileType::BINARY,
                             tk::real version = 2.2,
                             int datasize = sizeof(double) );

    //! Write Gmsh mesh to file
    void writeMesh( const UnsMesh& mesh );

  private:
    //! Write "$Nodes--$EndNodes" section
    void writeNodes( const UnsMesh& mesh );

    //! Write "$Elements--$EndElements" section
    void writeElements( const UnsMesh& mesh );

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

    //! Write element block: element ids and connectivity (node list)
    void writeElemBlock( std::size_t nnpe,
                         GmshElemType type,
                         const std::vector< std::size_t >& inpoel );

    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary
};

} // tk::

#endif // GmshMeshWriter_h
