// *****************************************************************************
/*!
  \file      src/IO/GmshMeshReader.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Gmsh mesh reader class declaration
  \details   Gmsh mesh reader class declaration. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
// *****************************************************************************
#ifndef GmshMeshReader_h
#define GmshMeshReader_h

#include <iosfwd>
#include <map>

#include "Types.h"
#include "Reader.h"
#include "GmshMeshIO.h"
#include "Exception.h"

namespace tk {

class UnsMesh;

//! Gmsh mesh reader
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   the Gmsh mesh generator: http://geuz.org/gmsh.
//! \author J. Bakosi
class GmshMeshReader : public Reader {

  public:
    //! Constructor
    explicit GmshMeshReader( const std::string filename ) :
      Reader(filename),
      m_version( 0.0 ),                        // 0.0: uninitialized
      m_datasize( 0 ),                         //   0: uninitialized
      m_type( GmshFileType::UNDEFINED )        //  -1: uninitialized
      {}

    //! Read Gmsh mesh
    void readMesh( UnsMesh& mesh );

  private:
    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes( UnsMesh& mesh );

    //! Read "$Elements--$EndElements" section
    void readElements( UnsMesh& mesh );

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames() __attribute__ ((noreturn));

    //! \brief Mesh ASCII type query
    //! \return true if member variable m_type indicates an ASCII mesh format
    //! \author J. Bakosi
    bool isASCII() const {
      Assert( m_type != GmshFileType::UNDEFINED, "Mesh type is undefined");
      return m_type == GmshFileType::ASCII ? true : false;
    }
    //! \brief Mesh binary type query
    //! \return true if member variable m_type indicates an binary mesh format
    //! \author J. Bakosi
    bool isBinary() const {
      Assert( m_type != GmshFileType::UNDEFINED, "Mesh type is undefined");
      return m_type == GmshFileType::BINARY ? true : false;
    }

    tk::real m_version;                 //!< Mesh version in mesh file
    int m_datasize;                     //!< Data size in mesh file
    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary

    //! \brief Gmsh element types and their corrseponding number of nodes
    //! \details See Gmsh documentation for element ids as keys
    const std::map< int, int > m_elemNodes {
      { GmshElemType::LIN,   2 },  // 2-node line
      { GmshElemType::TRI,   3 },  // 3-node triangle
      { GmshElemType::TET,   4 },  // 4-node tetrahedron
      { GmshElemType::PNT,   1 }   // 1-node point
    };
};

} // tk::

#endif // GmshMeshReader_h
