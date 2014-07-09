//******************************************************************************
/*!
  \file      src/IO/GmshMeshReader.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 08:58:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
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

//! GmshMeshReader : Reader
class GmshMeshReader : public tk::Reader {

  public:
    //! Constructor
    explicit GmshMeshReader( const std::string filename, UnsMesh& mesh ) :
      Reader(filename),
      m_version( 0.0 ),                        // 0.0: uninitialized
      m_datasize( 0 ),                         //   0: uninitialized
      m_type( GmshFileType::UNDEFINED ),       //  -1: uninitialized
      m_mesh( mesh ) {}

    //! Destructor, default compiler generated
    ~GmshMeshReader() noexcept override = default;

    //! Read Gmsh mesh
    void read() override;

  private:
    //! Don't permit copy constructor
    GmshMeshReader(const GmshMeshReader&) = delete;
    //! Don't permit copy assigment
    GmshMeshReader& operator=(const GmshMeshReader&) = delete;
    //! Don't permit move constructor
    GmshMeshReader(GmshMeshReader&&) = delete;
    //! Don't permit move assigment
    GmshMeshReader& operator=(GmshMeshReader&&) = delete;

    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();

    // Get mesh type queries
    bool isASCII() const {
      Assert( m_type != GmshFileType::UNDEFINED, "Mesh type is undefined");
      return m_type == GmshFileType::ASCII ? true : false;
    }
    bool isBinary() const {
      Assert( m_type != GmshFileType::UNDEFINED, "Mesh type is undefined");
      return m_type == GmshFileType::BINARY ? true : false;
    }

    tk::real m_version;                 //!< Mesh version in mesh file
    int m_datasize;                     //!< Data size in mesh file
    GmshFileType m_type;                //!< Mesh file type: 0:ASCII, 1:binary
    UnsMesh& m_mesh;                    //!< Mesh object

    //! Push back p for different element types
    template< class Pushed, class ElmType1, class ElmType2, class ElmType3 >
    void push_back( int etype, const Pushed& p, ElmType1& e1, ElmType2& e2,
                    ElmType3& e3 )
    {
      using tk::operator+;
      switch ( etype ) {
        case GmshElemType::LIN: e1.push_back( p ); break;
        case GmshElemType::TRI: e2.push_back( p ); break;
        case GmshElemType::TET: e3.push_back( p ); break;
        case GmshElemType::PNT: break;     // ignore 1-node 'point element' type
        default: Throw( std::string("Unsupported element type ") + etype +
                        " in mesh file: " + m_filename );

      }
    }

    //! Gmsh element types and their number of nodes,
    //! See Gmsh documentation for element ids as keys
    const std::map< int, int > m_elemNodes {
      { GmshElemType::LIN,   2 },  // 2-node line
      { GmshElemType::TRI,   3 },  // 3-node triangle
      { GmshElemType::TET,   4 },  // 4-node tetrahedron
      { GmshElemType::PNT,   1 }   // 1-node point
    };
};

} // quinoa::

#endif // GmshMeshReader_h
