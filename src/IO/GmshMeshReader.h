//******************************************************************************
/*!
  \file      src/IO/GmshMeshReader.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 08:26:20 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gmsh mesh reader class declaration
  \details   Gmsh mesh reader class declaration. Currently, this class supports
    line, triangle, tetrahedron, and point Gmsh element types.
*/
//******************************************************************************
#ifndef GmshMeshReader_h
#define GmshMeshReader_h

#include <map>

#include <Reader.h>
#include <UnsMesh.h>
#include <GmshMeshIO.h>
#include <Exception.h>

namespace quinoa {

//! \brief GmshMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a file saved by
//!   the Gmsh mesh generator: http://geuz.org/gmsh.
//! \author J. Bakosi
class GmshMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit GmshMeshReader( const std::string filename, UnsMesh& mesh ) :
      Reader(filename),
      m_version( 0.0 ),                        // 0.0: uninitialized
      m_datasize( 0 ),                         //   0: uninitialized
      m_type( GmshFileType::UNDEFINED ),       //  -1: uninitialized
      m_mesh( mesh ) {}

    //! Read Gmsh mesh
    void read() override;

  private:
    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

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
    UnsMesh& m_mesh;                    //!< Mesh object

    //! \brief Push back object p into underlying std::vector for different
    //!   element types
    template< class Pushed, class ElmType1, class ElmType2, class ElmType3 >
    void push_back( int etype,
                    const Pushed& p,
                    ElmType1& e1,
                    ElmType2& e2,
                    ElmType3& e3 )
    {
      using tk::operator<<;
      switch ( etype ) {
        case GmshElemType::LIN: e1.push_back( p ); break;
        case GmshElemType::TRI: e2.push_back( p ); break;
        case GmshElemType::TET: e3.push_back( p ); break;
        case GmshElemType::PNT: break;     // ignore 1-node 'point element' type
        default: Throw( std::string("Unsupported element type ") << etype <<
                        " in mesh file: " << m_filename );

      }
    }

    //! \brief Gmsh element types and their corrseponding number of nodes
    //! \details See Gmsh documentation for element ids as keys
    const std::map< int, int > m_elemNodes {
      { GmshElemType::LIN,   2 },  // 2-node line
      { GmshElemType::TRI,   3 },  // 3-node triangle
      { GmshElemType::TET,   4 },  // 4-node tetrahedron
      { GmshElemType::PNT,   1 }   // 1-node point
    };
};

} // quinoa::

#endif // GmshMeshReader_h
